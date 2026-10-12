"""Crash-safe checkpoints; exclusive creation protects previous run receipts."""
import json
import os
from pathlib import Path
import tempfile


def reserve(path, receipt):
    content = json.dumps(receipt, indent=2, default=str) + '\n'
    with Path(path).open('x') as stream:
        stream.write(content)


def save(path, receipt):
    path = Path(path)
    content = json.dumps(receipt, indent=2, default=str) + '\n'
    fd, name = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_partial(path):
    """Provider wrappers must preserve timeout status even with no valid child file."""
    try:
        value = json.loads(Path(path).read_text())
        if not isinstance(value, dict):
            raise ValueError('Receipt root must be an object')
        return value
    except (OSError, ValueError) as error:
        return {'status': 'receipt_unavailable', 'gpu_executed': None,
                'error_type': type(error).__name__, 'error': str(error)}
