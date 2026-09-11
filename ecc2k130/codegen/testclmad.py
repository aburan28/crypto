"""Check CLMAD selection and compiler guards without requiring a CUDA installation."""
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ClmadGuardTests(unittest.TestCase):
    def preprocess(self, *flags):
        compiler = shlex.split(os.environ.get('CXX', 'c++'))
        if not compiler or shutil.which(compiler[0]) is None:
            self.skipTest('C++ preprocessor unavailable')
        return subprocess.run(
            compiler + ['-std=c++17', '-E', '-P', '-x', 'c++', '-I', str(ROOT / 'include'),
                        *flags, '-'], input='#include "packed131.h"\n',
            text=True, capture_output=True, timeout=30)

    def test_default_and_host_fallback_do_not_emit_native_instructions(self):
        default = self.preprocess()
        disabled = self.preprocess('-DECC_PACKED_CLMAD=0')
        host = self.preprocess('-DECC_PACKED_CLMAD=1')
        for result in (default, disabled, host):
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn('clmad.lo.u64', result.stdout)
        self.assertEqual(default.stdout, disabled.stdout)

    def test_invalid_flag_is_rejected(self):
        for value in ('-1', '2'):
            result = self.preprocess('-DECC_PACKED_CLMAD=' + value)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('ECC_PACKED_CLMAD must be 0 or 1', result.stderr)

    def test_unsupported_compiler_or_device_is_rejected(self):
        # These are preprocessing controls, not claims of CUDA compilation.
        common = ['-DECC_PACKED_CLMAD=1', '-D__CUDACC__', '-D__CUDACC_VER_MAJOR__=13']
        old_compiler = self.preprocess(*common, '-D__CUDACC_VER_MINOR__=2', '-D__CUDA_ARCH__=800')
        old_device = self.preprocess(*common, '-D__CUDACC_VER_MINOR__=3', '-D__CUDA_ARCH__=750')
        self.assertNotEqual(old_compiler.returncode, 0)
        self.assertIn('requires CUDA 13.3', old_compiler.stderr)
        self.assertNotEqual(old_device.returncode, 0)
        self.assertIn('requires sm_80', old_device.stderr)

    def test_supported_preprocessor_configuration_selects_both_halves(self):
        result = self.preprocess('-DECC_PACKED_CLMAD=1', '-D__CUDACC__',
                                 '-D__CUDACC_VER_MAJOR__=13', '-D__CUDACC_VER_MINOR__=3',
                                 '-D__CUDA_ARCH__=800')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('clmad.lo.u64', result.stdout)
        self.assertIn('clmad.hi.u64', result.stdout)


if __name__ == '__main__':
    unittest.main()
