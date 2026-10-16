"""No-signal guards for the bounded Docker recovery attempt."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

import recover


class GracefulRecoveryGuardTest(unittest.TestCase):
    def test_changed_pid_identity_blocks_signal_and_output(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "new-output"
            wrong = {"pid": recover.PARENT_PID, "ppid": 1,
                     "state": "S", "command": "/usr/bin/other"}
            with mock.patch.object(recover, "process", return_value=wrong), \
                 mock.patch.object(recover, "ping", return_value=(
                     "HTTP/1.1 503 Service Unavailable\r\n\r\n"
                     "Docker Desktop is unable to start")), \
                 mock.patch.object(recover.os, "kill") as kill:
                with self.assertRaisesRegex(ValueError, "signal withheld"):
                    recover.run(output)
                kill.assert_not_called()
                self.assertFalse(output.exists())

    def test_healthy_socket_blocks_signal_and_output(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "new-output"
            parent = {"pid": recover.PARENT_PID, "ppid": 1,
                      "state": "S", "command": recover.PARENT_COMMAND}
            with mock.patch.object(recover, "process", return_value=parent), \
                 mock.patch.object(recover, "ping", return_value="HTTP/1.1 200 OK"), \
                 mock.patch.object(recover.os, "kill") as kill:
                with self.assertRaisesRegex(ValueError, "signal withheld"):
                    recover.run(output)
                kill.assert_not_called()
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
