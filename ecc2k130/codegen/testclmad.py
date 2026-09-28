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
        # Old toolkits still hard-error. sm_75 in a fat ARCHES="75 …" CLMAD=1
        # build must soft-fall back to the software product so g4dn stays
        # buildable (packed131.h ECC_USE_CLMAD_INSN).
        common = ['-DECC_PACKED_CLMAD=1', '-D__CUDACC__', '-D__CUDACC_VER_MAJOR__=13']
        old_compiler = self.preprocess(*common, '-D__CUDACC_VER_MINOR__=2', '-D__CUDA_ARCH__=800')
        old_device = self.preprocess(*common, '-D__CUDACC_VER_MINOR__=3', '-D__CUDA_ARCH__=750')
        self.assertNotEqual(old_compiler.returncode, 0)
        self.assertIn('requires CUDA 13.3', old_compiler.stderr)
        self.assertEqual(old_device.returncode, 0, old_device.stderr)
        self.assertNotIn('clmad.lo.u64', old_device.stdout)
        self.assertNotIn('clmad.hi.u64', old_device.stdout)

    def test_supported_preprocessor_configuration_selects_both_halves(self):
        result = self.preprocess('-DECC_PACKED_CLMAD=1', '-D__CUDACC__',
                                 '-D__CUDACC_VER_MAJOR__=13', '-D__CUDACC_VER_MINOR__=3',
                                 '-D__CUDA_ARCH__=800')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('clmad.lo.u64', result.stdout)
        self.assertIn('clmad.hi.u64', result.stdout)

    def test_square_selection_and_host_fallback(self):
        # Native squaring follows ECC_USE_CLMAD_INSN (CLMAD + sm_80+), not the
        # historical ECC_PACKED_CLMAD_SQUARE selector. Host and Turing keep the
        # shift/mask path; the square flag still validates 0/1 and CLMAD dependency.
        host = self.preprocess('-DECC_PACKED_CLMAD=1', '-DECC_PACKED_CLMAD_SQUARE=0')
        self.assertEqual(host.returncode, 0, host.stderr)
        self.assertNotIn('clmad.lo.u64', host.stdout)
        device = self.preprocess('-DECC_PACKED_CLMAD=1', '-DECC_PACKED_CLMAD_SQUARE=0',
                                 '-D__CUDACC__', '-D__CUDACC_VER_MAJOR__=13',
                                 '-D__CUDACC_VER_MINOR__=3', '-D__CUDA_ARCH__=1200')
        self.assertEqual(device.returncode, 0, device.stderr)
        self.assertIn('clmad.lo.u64 %0, %1, %1, 0;', device.stdout)
        for flags in (('-DECC_PACKED_CLMAD_SQUARE=1',),
                      ('-DECC_PACKED_CLMAD=1', '-DECC_PACKED_CLMAD_SQUARE=2')):
            result = self.preprocess(*flags)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('ECC_PACKED_CLMAD_SQUARE', result.stderr)


class NativeSquareProofTests(unittest.TestCase):
    def test_exact_source_equivalence(self):
        from prove_native_square import prove, HEADER
        result = prove(HEADER.read_text())
        self.assertTrue(result['equivalent'])
        self.assertEqual(result['covered_inputs'], 2**32)
        self.assertEqual(result['nonzero_output_bits'], 32)

    def test_wrong_mask_and_changed_native_operands_are_detected(self):
        from prove_native_square import prove, HEADER
        source = HEADER.read_text()
        bad_mask = source.replace('0x5555555555555555ull;', '0x5555555555555554ull;')
        with self.assertRaisesRegex(AssertionError, 'inequivalent output bits'):
            prove(bad_mask)
        for before, after in [('clmad.lo.u64 %0, %1, %1, 0;', 'clmad.hi.u64 %0, %1, %1, 0;'),
                              ('"l"((uint64_t)x)', '"l"((uint64_t)(x >> 1))')]:
            with self.assertRaisesRegex(ValueError, 'native branch'):
                prove(source.replace(before, after))


if __name__ == '__main__':
    unittest.main()

