"""Exercise decoded equivalence, worker independence and benchmark failures."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = sys.argv.pop(1)
LAUNCHER = str(Path(BINARY).with_name('xml_benchmark_launcher'))


def invoke(*args):
    result = subprocess.run([LAUNCHER, BINARY, *map(str, args)], check=True, capture_output=True, text=True)
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
                for mode in ('stream', 'load', 'indexed', 'index-open', 'contexts', 'idle', 'indexed-metadata', 'contexts-metadata'):
                    single = invoke(mode, path, 1, 41)
                    parallel = invoke(mode, path, 4, 41)
                    self.assertEqual(single['hash_per_worker'], parallel['hash_per_worker'])
                    self.assertEqual(single['points_per_worker'], parallel['points_per_worker'])
                    if mode in ('stream', 'load'):
                        self.assertEqual(single['hash_per_worker'], result['hash'])
                    if mode.startswith('contexts'):
                        self.assertEqual(single['index_spectra'], 17)
                        self.assertGreater(single['resident_rss_bytes'], 0)
                        self.assertGreater(parallel['resident_rss_bytes'], 0)
                    if mode == 'indexed-metadata':
                        self.assertEqual(single['hash_per_worker'], invoke('indexed', path, 1, 41)['hash_per_worker'])

    def test_errors_propagate(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'broken.mzML'
            path.write_text('<mzML><broken></mzML>')
            for mode in ('stream', 'load', 'indexed', 'contexts', 'indexed-metadata', 'contexts-metadata'):
                result = subprocess.run([LAUNCHER, BINARY, mode, str(path), '2', '3'], capture_output=True)
                self.assertNotEqual(result.returncode, 0)

    def test_stale_registration_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'registration.mzML'
            invoke('generate', path, 1, 3, 'plain')
            code = '''import os, sys
p = '/tmp/__KMP_REGISTERED_LIB_' + str(os.getpid())
try:
    fd = os.open(p, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
except FileExistsError:
    sys.exit(77)
os.close(fd)
os.utime(p, (0, 0))
os.execv(sys.argv[1], sys.argv[1:])
'''
            result = subprocess.run([sys.executable, '-c', code, LAUNCHER, BINARY,
                                     'stream', str(path), '1', '1'], capture_output=True, text=True)
            if result.returncode == 77:
                self.skipTest('Own PID already has an existing runtime registration file')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('XML_BENCHMARK_REMOVED_STALE_KMP', result.stderr)


if __name__ == '__main__':
    unittest.main()
