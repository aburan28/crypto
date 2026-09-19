#!/usr/bin/env python3
"""Operator entry point: `python3 -m controlplane <command>`.

    migrate     apply the schema (idempotent; safe on every boot)
    status      one JSON object: slots, corpus, fleet, cache
    slots       the registry, or one slot's event history
    doctor      check each store and say which one is unhappy
    prune       drop worker rows and slot events older than a cutoff
    retire      take a run id out of service

Every command prints one JSON object, like `ca` and like `worker.py`'s status
output, so the dashboard, the alarms and a person all read the same thing.
`doctor` is the one to run first when something is wrong: it reports each of
the three stores separately, because "the campaign is down" is almost always
one of them being down and the other two being fine.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

if __package__ in (None, ""):      # running the file directly
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from controlplane.cache import cacheFromEnv
    from controlplane.config import Config
    from controlplane.coordinator import Coordinator
    from controlplane.db import databaseFromEnv
    from controlplane.schema import migrate, pending
    from controlplane.store import storeFromEnv
else:
    from .cache import cacheFromEnv
    from .config import Config
    from .coordinator import Coordinator
    from .db import databaseFromEnv
    from .schema import migrate, pending
    from .store import storeFromEnv


def emit(obj):
    json.dump(obj, sys.stdout, indent=2, sort_keys=True, default=str)
    sys.stdout.write("\n")


def doctor(config):
    """Report each store separately, and never let one failure hide another."""
    out = {"campaign": config.campaign, "at": int(time.time())}

    out["database"] = {"configured": bool(config.databaseUrl or config.dbHost or config.localStore)}
    try:
        db = databaseFromEnv(config)
        db.fetchOne("SELECT 1")
        out["database"].update(ok=True, pending=pending(db), dialect=db.dialect)
    except Exception as exc:
        out["database"].update(ok=False, error="%s: %s" % (type(exc).__name__, exc))

    out["cache"] = {"configured": bool(config.redisUrl) and config.cacheEnabled}
    try:
        cache = cacheFromEnv(config)
        if not cache.enabled:
            # Not an error: a fleet with no cache is a supported configuration
            # and is the control any cache claim is measured against.
            out["cache"].update(ok=True, note="no cluster configured; RDS answers everything")
        else:
            cache.setJson("doctor", {"at": int(time.time())}, ttl=60)
            out["cache"].update(ok=cache.getJson("doctor") is not None, stats=cache.stats())
    except Exception as exc:
        out["cache"].update(ok=False, error="%s: %s" % (type(exc).__name__, exc))

    out["corpus"] = {"configured": bool(config.bucket or config.localStore)}
    try:
        store = storeFromEnv(config)
        store.list(limit=1)
        out["corpus"].update(ok=True, bucket=getattr(store, "bucket", getattr(store, "root", "")))
    except Exception as exc:
        out["corpus"].update(ok=False, error="%s: %s" % (type(exc).__name__, exc))

    out["ok"] = all(section.get("ok") for section in
                    (out["database"], out["cache"], out["corpus"]))
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--campaign", default=None, help="override RHO_CAMPAIGN")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("migrate", help="apply the control-plane schema")
    sub.add_parser("status", help="slots, corpus, fleet and cache in one object")
    sub.add_parser("doctor", help="check RDS, ElastiCache and S3 separately")

    p = sub.add_parser("slots", help="the slot registry")
    p.add_argument("--slot", type=int, default=None)
    p.add_argument("--events", action="store_true", help="show the event log instead")
    p.add_argument("--limit", type=int, default=100)

    p = sub.add_parser("fleet", help="live workers")
    p.add_argument("--fresh", action="store_true", help="skip the cache")

    p = sub.add_parser("prune", help="drop stale worker rows and slot events")
    p.add_argument("--days", type=float, default=7.0)

    p = sub.add_parser("retire", help="take a run id out of service")
    p.add_argument("slot", type=int)
    p.add_argument("--reason", required=True)
    p.add_argument("--state", default="retired", choices=("retired", "solved", "error"))

    args = parser.parse_args(argv)
    config = Config.fromEnv()
    if args.campaign:
        config = Config(**dict(config.__dict__, campaign=args.campaign))

    if args.command == "doctor":
        out = doctor(config)
        emit(out)
        return 0 if out["ok"] else 1

    if args.command == "migrate":
        db = databaseFromEnv(config)
        emit({"applied": migrate(db), "pending": pending(db)})
        return 0

    coordinator = Coordinator.fromEnv(config)
    try:
        if args.command == "status":
            emit(coordinator.status())
        elif args.command == "slots":
            if args.events:
                emit({"events": coordinator.slots.events(args.slot, limit=args.limit)})
            elif args.slot is not None:
                emit(coordinator.slots.get(args.slot) or {})
            else:
                emit({"slots": coordinator.slots.scan(fresh=True)})
        elif args.command == "fleet":
            emit(coordinator.fleet(fresh=args.fresh))
        elif args.command == "prune":
            emit(coordinator.prune(before=int(time.time() - args.days * 86400)))
        elif args.command == "retire":
            ok = coordinator.slots.retire(args.slot, args.reason, state=args.state)
            emit({"slot": args.slot, "state": args.state, "changed": ok})
            return 0 if ok else 1
    finally:
        coordinator.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
