import os
import sys

# Put the tpu/ directory on the path so `import ic.*` resolves the package.
_TPU_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _TPU_DIR not in sys.path:
    sys.path.insert(0, _TPU_DIR)

# JAX on CPU; keep it single-threaded and quiet for deterministic tests.
os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=1")
