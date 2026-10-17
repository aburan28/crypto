# Darwin memory-limit attempts

Before the 60-child V1/V2 panels, the native process meter attempted
to set a 536,870,912-byte `RLIMIT_AS` on the solver child. The
`Command::status().unwrap()` call failed while spawning the child:
`Os { code: 22, kind: InvalidInput, message: "Invalid argument" }`.
The meter process exited 101. The empty child stdout/stderr files
were preserved as `METER_FAILED_RLIMIT_AS_V1.*`.

A second version tried the same value under `RLIMIT_DATA` on Darwin
and failed at `Command::status().unwrap()` with the same code 22
`Invalid argument`; it also exited 101. Its empty child files are
`METER_FAILED_RLIMIT_DATA_V1.*`. These failures concern host
resource-limit enforcement, not χ decomposition or group arithmetic.

The subsequently used meter enforces a 300-second wall limit, records
whole-child user/system CPU and observed peak RSS, and fails a child
receipt if observed RSS exceeds 536,870,912 bytes. It reports
`rlimit_as_applied=false` and `rlimit_data_applied=false` on
Darwin. On Linux it applies `RLIMIT_AS` before child execution.
