//! Small, streaming SQLite adapter for the Linux checkpoint example.
//! All handles and column pointers stay inside the owning connection lifetime.
use std::ffi::{c_char, c_int, c_void, CStr, CString};
use std::marker::PhantomData;

#[link(name = "libsqlite3.so.0", kind = "dylib", modifiers = "+verbatim")]
unsafe extern "C" {
    fn sqlite3_open(filename: *const c_char, db: *mut *mut c_void) -> c_int;
    fn sqlite3_close(db: *mut c_void) -> c_int;
    fn sqlite3_errmsg(db: *mut c_void) -> *const c_char;
    fn sqlite3_exec(
        db: *mut c_void,
        sql: *const c_char,
        cb: *const c_void,
        arg: *mut c_void,
        err: *mut *mut c_char,
    ) -> c_int;
    fn sqlite3_prepare_v2(
        db: *mut c_void,
        sql: *const c_char,
        n: c_int,
        stmt: *mut *mut c_void,
        tail: *mut *const c_char,
    ) -> c_int;
    fn sqlite3_step(stmt: *mut c_void) -> c_int;
    fn sqlite3_finalize(stmt: *mut c_void) -> c_int;
    fn sqlite3_column_int64(stmt: *mut c_void, col: c_int) -> i64;
    fn sqlite3_column_text(stmt: *mut c_void, col: c_int) -> *const u8;
    fn sqlite3_busy_timeout(db: *mut c_void, ms: c_int) -> c_int;
    fn sqlite3_db_status(
        db: *mut c_void,
        op: c_int,
        current: *mut c_int,
        highwater: *mut c_int,
        reset: c_int,
    ) -> c_int;
}

pub fn quote(s: &str) -> String {
    format!("'{}'", s.replace('\'', "''"))
}

pub struct Db {
    raw: *mut c_void,
}
pub struct Statement<'a> {
    raw: *mut c_void,
    db: &'a Db,
    _not_send: PhantomData<*mut ()>,
}

impl Db {
    pub fn open(path: &str) -> Self {
        let name = CString::new(path).expect("database path contains NUL");
        let mut raw = std::ptr::null_mut();
        let code = unsafe { sqlite3_open(name.as_ptr(), &mut raw) };
        if code != 0 {
            let msg = if raw.is_null() {
                "no SQLite handle".to_owned()
            } else {
                unsafe { CStr::from_ptr(sqlite3_errmsg(raw)) }
                    .to_string_lossy()
                    .into_owned()
            };
            if !raw.is_null() {
                unsafe {
                    sqlite3_close(raw);
                }
            }
            panic!("SQLite open: {msg}");
        }
        let db = Self { raw };
        db.check(unsafe { sqlite3_busy_timeout(raw, 10_000) });
        db.exec("PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL; PRAGMA cache_size=-8192; PRAGMA temp_store=FILE; PRAGMA mmap_size=0; PRAGMA foreign_keys=ON;");
        db
    }
    fn check(&self, code: c_int) {
        assert_eq!(
            code,
            0,
            "SQLite: {}",
            unsafe { CStr::from_ptr(sqlite3_errmsg(self.raw)) }.to_string_lossy()
        );
    }
    pub fn exec(&self, sql: &str) {
        let sql = CString::new(sql).expect("SQL contains NUL");
        self.check(unsafe {
            sqlite3_exec(
                self.raw,
                sql.as_ptr(),
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            )
        });
    }
    pub fn prepare(&self, sql: &str) -> Statement<'_> {
        let sql = CString::new(sql).expect("SQL contains NUL");
        let mut raw = std::ptr::null_mut();
        self.check(unsafe {
            sqlite3_prepare_v2(self.raw, sql.as_ptr(), -1, &mut raw, std::ptr::null_mut())
        });
        assert!(!raw.is_null());
        Statement {
            raw,
            db: self,
            _not_send: PhantomData,
        }
    }
    pub fn int(&self, sql: &str) -> i64 {
        let mut s = self.prepare(sql);
        assert!(s.step(), "missing scalar row");
        s.int(0)
    }
    pub fn text(&self, sql: &str) -> String {
        let mut s = self.prepare(sql);
        assert!(s.step(), "missing scalar row");
        s.text(0)
    }
    pub fn cache_bytes(&self) -> i32 {
        let (mut current, mut peak) = (0, 0);
        self.check(unsafe { sqlite3_db_status(self.raw, 1, &mut current, &mut peak, 0) });
        current
    }
}
impl Drop for Db {
    fn drop(&mut self) {
        unsafe {
            sqlite3_close(self.raw);
        }
    }
}
impl Statement<'_> {
    pub fn step(&mut self) -> bool {
        match unsafe { sqlite3_step(self.raw) } {
            100 => true,
            101 => false,
            code => {
                self.db.check(code);
                unreachable!()
            }
        }
    }
    pub fn int(&self, col: i32) -> i64 {
        unsafe { sqlite3_column_int64(self.raw, col) }
    }
    pub fn text(&self, col: i32) -> String {
        let p = unsafe { sqlite3_column_text(self.raw, col) };
        assert!(!p.is_null(), "unexpected NULL text");
        unsafe { CStr::from_ptr(p.cast()) }
            .to_str()
            .expect("UTF-8 column")
            .to_owned()
    }
}
impl Drop for Statement<'_> {
    fn drop(&mut self) {
        unsafe {
            sqlite3_finalize(self.raw);
        }
    }
}
