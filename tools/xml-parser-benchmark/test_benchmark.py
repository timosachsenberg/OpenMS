"""Exercise decoded equivalence, worker independence and benchmark failures."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = sys.argv.pop(1)


def invoke(*args):
    result = subprocess.run([BINARY, *map(str, args)], check=True, capture_output=True, text=True)
    lines = [line for line in result.stdout.splitlines() if line.startswith('{')]
    return json.loads(lines[-1]) if lines else None


class BenchmarkTest(unittest.TestCase):
    def test_equivalence_and_workers(self):
        with tempfile.TemporaryDirectory() as directory:
            for compression in ('plain', 'zlib'):
                path = Path(directory) / f'{compression}.mzML'
                invoke('generate', path, 17, 31, compression)
                result = invoke('verify', path)
                self.assertTrue(result['verified'])
                self.assertEqual(result['spectra'], 17)
                self.assertEqual(result['points'], 17 * 31 + 17)
                for mode in ('stream', 'indexed', 'index-open'):
                    single = invoke(mode, path, 1, 41)
                    parallel = invoke(mode, path, 4, 41)
                    self.assertEqual(single['hash_per_worker'], parallel['hash_per_worker'])
                    self.assertEqual(single['points_per_worker'], parallel['points_per_worker'])

    def test_errors_propagate(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'broken.mzML'
            path.write_text('<mzML><broken></mzML>')
            for mode in ('stream', 'indexed'):
                result = subprocess.run([BINARY, mode, str(path), '2', '3'], capture_output=True)
                self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
