"""RDS Postgres access: one connection per thread, retried across a failover.

A Multi-AZ failover is a normal event, not an outage: the writer endpoint
resolves to the new instance within about a minute and every connection open
across the switch is dead.  A control plane that treats that as a crash takes
the whole fleet down with it, so every statement here runs through `run()`,
which reconnects and retries on a connection-shaped error and does not retry
anything else.  A syntax error or a constraint violation is a bug and is
raised immediately.

Retries are safe because every statement this package issues is idempotent or
conditional: claims and heartbeats are conditional UPDATEs that a replay
either wins once or loses, ledger inserts are `ON CONFLICT DO NOTHING`, and
nothing accumulates.  A caller adding a bare `UPDATE ... SET x = x + 1` would
break that and must wrap it in `transaction()` instead, which does not retry.

`psycopg` (v3) is imported lazily, inside `connect()`.  Tests drive this class
with a `connect` callable of their own -- SQLite, in the suite -- so the SQL
in `schema.py` and `slots.py` is exercised for real without a database.
"""

from __future__ import annotations

import json
import os
import random
import threading
import time
import urllib.parse

from .config import Config


class DatabaseError(RuntimeError):
    pass


# psycopg raises these for a connection that died under a statement; the
# classes are matched by name so this module does not import psycopg to
# decide whether to retry.  SQLState 57P01 (admin shutdown) and 08xxx
# (connection exception) are the failover's own codes.
_RETRY_CLASSES = ("OperationalError", "InterfaceError", "AdminShutdown",
                  "ConnectionException", "ConnectionDoesNotExist",
                  "ConnectionFailure", "CannotConnectNow")
_RETRY_SQLSTATES = ("57P01", "57P02", "57P03", "40001", "40P01")


def isTransient(exc):
    """True for an error a reconnect could plausibly fix.

    Name-based rather than type-based on purpose: this module must import
    without psycopg, and the suite drives it with a different driver whose
    exception hierarchy is unrelated.
    """
    if isinstance(exc, DatabaseError):
        return False
    state = getattr(exc, "sqlstate", None) or getattr(exc, "pgcode", None)
    if state and str(state) in _RETRY_SQLSTATES:
        return True
    if state:
        return str(state).startswith("08")
    names = {c.__name__ for c in type(exc).__mro__}
    if names & set(_RETRY_CLASSES):
        return True
    text = str(exc).lower()
    return any(s in text for s in (
        "server closed the connection", "connection already closed",
        "terminating connection", "could not connect", "connection reset",
        "ssl connection has been closed", "consuming input failed"))


def databaseUrl(config=None):
    """The libpq URL, from the environment or from Secrets Manager.

    Identical resolution order to `dp_ingest.py`, deliberately: the ingest and
    the control plane must never end up pointed at two different databases
    because one of them preferred a different variable.
    """
    config = config or Config.fromEnv()
    if config.databaseUrl:
        return config.databaseUrl
    if not config.dbHost:
        raise DatabaseError("set DATABASE_URL, or RHO_DB_HOST plus RHO_DB_SECRET")
    import boto3  # lazy: this module imports on a host with no AWS SDK

    blob = boto3.client("secretsmanager").get_secret_value(SecretId=config.dbSecret)["SecretString"]
    c = json.loads(blob)
    name = config.dbName or c.get("dbname") or "postgres"
    port = c.get("port") or 5432
    return "postgresql://%s:%s@%s:%s/%s?sslmode=%s" % (
        urllib.parse.quote(c["username"], safe=""),
        urllib.parse.quote(c["password"], safe=""),
        config.dbHost, port, name, config.dbSslMode)


class Database:
    """A thread-safe handle on one Postgres database.

    Connections are per thread and never shared: psycopg connections are not
    thread-safe, and a shared one under the fleet's heartbeat rate is a
    deadlock waiting for a quiet afternoon.
    """

    def __init__(self, connect, paramstyle="format", dialect="postgres",
                 attempts=4, sleep=time.sleep):
        self._connect = connect
        self.paramstyle = paramstyle      # "format" (%s) or "qmark" (?)
        self.dialect = dialect            # "postgres" or "sqlite"
        self.attempts = max(1, attempts)
        self._sleep = sleep
        self._local = threading.local()
        self._lock = threading.Lock()
        self._connections = []            # every connection handed out, for close()

    # -- construction -------------------------------------------------------
    @classmethod
    def fromEnv(cls, config=None):
        config = config or Config.fromEnv()
        url = databaseUrl(config)

        def connect():
            import psycopg  # lazy

            conn = psycopg.connect(url, connect_timeout=config.dbConnectTimeout,
                                   autocommit=False)
            # A control-plane statement that runs longer than this is stuck,
            # and a stuck claim holds a worker's GPU idle.  Fail instead.
            with conn.cursor() as cur:
                cur.execute("SET statement_timeout = %s" % int(config.dbStatementTimeoutMs))
                cur.execute("SET idle_in_transaction_session_timeout = %s"
                            % int(max(config.dbStatementTimeoutMs, 30000)))
            conn.commit()
            return conn

        return cls(connect)

    # -- connection handling ------------------------------------------------
    def connection(self):
        conn = getattr(self._local, "conn", None)
        if conn is None:
            conn = self._connect()
            self._local.conn = conn
            with self._lock:
                self._connections.append(conn)
        return conn

    def _drop(self):
        conn = getattr(self._local, "conn", None)
        self._local.conn = None
        if conn is None:
            return
        with self._lock:
            if conn in self._connections:
                self._connections.remove(conn)
        try:
            conn.close()
        except Exception:
            pass

    def close(self):
        self._local.conn = None
        with self._lock:
            conns, self._connections = self._connections, []
        for conn in conns:
            try:
                conn.close()
            except Exception:
                pass

    # -- statements ---------------------------------------------------------
    def sql(self, statement):
        """Translate the package's `?` placeholders for the driver in use."""
        if self.paramstyle == "qmark":
            return statement
        return statement.replace("?", "%s")

    def _once(self, statements, fetch):
        conn = self.connection()
        rows = None
        try:
            # Not `with conn.cursor() as cur`: psycopg's cursor is a context
            # manager and sqlite3's is not, and this module runs on both.
            cur = conn.cursor()
            try:
                for stmt, params in statements:
                    cur.execute(self.sql(stmt), tuple(params))
                    if fetch == "all":
                        rows = cur.fetchall()
                    elif fetch == "one":
                        rows = cur.fetchone()
                    elif fetch == "count":
                        rows = cur.rowcount
            finally:
                cur.close()
            conn.commit()
            return rows
        except Exception:
            try:
                conn.rollback()
            except Exception:
                self._drop()
            raise

    def run(self, statement, params=(), fetch=None):
        """One statement, committed, retried across a failover.

        `fetch` is None, "one", "all" or "count" (the driver's rowcount, which
        is how a conditional UPDATE reports whether it won).
        """
        return self.runAll([(statement, params)], fetch=fetch)

    def runAll(self, statements, fetch=None):
        """Several statements in one transaction, retried as a unit."""
        last = None
        for attempt in range(self.attempts):
            try:
                return self._once(statements, fetch)
            except Exception as exc:
                last = exc
                if not isTransient(exc) or attempt == self.attempts - 1:
                    raise
                self._drop()
                # Jittered backoff: a failover wakes every worker in the
                # fleet at once, and an unjittered retry storm is how a
                # database that just came back goes away again.
                delay = min(8.0, 0.25 * (2 ** attempt))
                self._sleep(delay * (0.5 + random.random()))
        raise last

    def transaction(self):
        """A context manager over one connection's cursor, no retry.

        For the multi-statement sequences whose correctness depends on being
        one transaction -- the fenced writes in `coordinator.py`.  A caller
        that wants a retry on top of this must be able to replay the whole
        block, and must say so by calling it again itself.
        """
        return _Transaction(self)

    # -- conveniences -------------------------------------------------------
    def fetchAll(self, statement, params=()):
        return self.run(statement, params, fetch="all") or []

    def fetchOne(self, statement, params=()):
        return self.run(statement, params, fetch="one")

    def execute(self, statement, params=()):
        return self.run(statement, params, fetch="count")


class _Transaction:
    def __init__(self, db):
        self.db = db
        self.conn = None
        self.cur = None

    def __enter__(self):
        self.conn = self.db.connection()
        self.cur = self.conn.cursor()
        return self

    def execute(self, statement, params=()):
        self.cur.execute(self.db.sql(statement), tuple(params))
        return self.cur

    def fetchOne(self, statement, params=()):
        return self.execute(statement, params).fetchone()

    def fetchAll(self, statement, params=()):
        return self.execute(statement, params).fetchall()

    def __exit__(self, kind, value, tb):
        try:
            self.cur.close()
        except Exception:
            pass
        if kind is None:
            self.conn.commit()
        else:
            try:
                self.conn.rollback()
            except Exception:
                self.db._drop()
        return False


def sqliteDatabase(path=":memory:", **kw):
    """A Database over SQLite, for rehearsals and for the test suite.

    Every statement in this package is written to run on both engines: `?`
    placeholders, integer epoch timestamps rather than `now()`, no server-side
    sequences, and `ON CONFLICT` clauses SQLite understands.  That is not a
    portability goal for its own sake -- it is what lets the claim/fence logic
    be tested against a real SQL engine on a laptop with no database.
    """
    import sqlite3

    def connect():
        conn = sqlite3.connect(path, timeout=30, isolation_level="DEFERRED",
                               check_same_thread=False)
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA busy_timeout = 30000")
        conn.commit()
        return conn

    return Database(connect, paramstyle="qmark", dialect="sqlite", **kw)


def databaseFromEnv(config=None):
    """The database a CLI should use: SQLite under ECC_LOCAL_STORE, else RDS."""
    config = config or Config.fromEnv()
    if config.localStore:
        os.makedirs(config.localStore, exist_ok=True)
        return sqliteDatabase(os.path.join(config.localStore, "controlplane.db"))
    return Database.fromEnv(config)
