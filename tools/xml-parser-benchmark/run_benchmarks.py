"""Run identical OpenMS workloads in fresh processes and record time/RSS/digests.

Example: python run_benchmarks.py --binary build/xml_parser_benchmark \
  --backend xerces=/snapshots/xerces/lib --backend libxml2=/snapshots/libxml2/lib \
  --input data/medium.mzML --input data/large.mzML --output results.json
"""
import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time


def measured(command, environment, timeout):
    with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
        begin = time.perf_counter()
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=environment)
        deadline = begin + timeout
        while True:
            child, status, usage = os.wait4(process.pid, os.WNOHANG)
            if child:
                break
            if time.perf_counter() >= deadline:
                process.kill()
                child, status, usage = os.wait4(process.pid, 0)
                process.returncode = os.waitstatus_to_exitcode(status)
                raise TimeoutError(f'Benchmark timed out: {command}')
            time.sleep(0.02)
        process.returncode = os.waitstatus_to_exitcode(status)
        stdout.seek(0)
        stderr.seek(0)
        output, errors = stdout.read().decode(), stderr.read().decode()
        if process.returncode:
            raise RuntimeError(f'{command} exited {process.returncode}\n{output}\n{errors}')
        lines = [line for line in output.splitlines() if line.startswith('{')]
        result = json.loads(lines[-1])
        result.update(process_seconds=time.perf_counter() - begin,
                      max_rss_bytes=usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
                      user_seconds=usage.ru_utime, system_seconds=usage.ru_stime)
        return result


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def summarize(records):
    groups = defaultdict(list)
    for record in records:
        groups[(record['backend'], record['input'], record['mode'], record['workers'])].append(record)
    result = []
    for (backend, filename, mode, workers), rows in sorted(groups.items()):
        seconds = [row['seconds'] for row in rows]
        rss = [row['max_rss_bytes'] / 1024 ** 2 for row in rows]
        median = statistics.median(seconds)
        item = dict(backend=backend, input=filename, mode=mode, workers=workers,
                    repetitions=len(rows), seconds_median=median, seconds_min=min(seconds),
                    seconds_max=max(seconds), rss_mib_median=statistics.median(rss),
                    rss_mib_max=max(rss), hash_per_worker=rows[0]['hash_per_worker'])
        item['job_open_seconds_median'] = statistics.median(job['open_seconds'] for row in rows for job in row['jobs'])
        item['job_read_seconds_median'] = statistics.median(job['read_seconds'] for row in rows for job in row['jobs'])
        if mode == 'stream':
            item['aggregate_mib_per_second'] = rows[0]['file_bytes'] * workers / median / 1024 ** 2
        elif mode == 'indexed':
            item['aggregate_spectra_per_second'] = rows[0]['spectra_per_worker'] * workers / median
        result.append(item)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--backend', required=True, action='append', help='LABEL=directory of OpenMS libraries')
    parser.add_argument('--input', required=True, action='append', type=Path)
    parser.add_argument('--workers', default='1,2,4,8')
    parser.add_argument('--modes', default='stream,indexed,index-open')
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--queries', type=int, default=2048)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--output', required=True, type=Path)
    arguments = parser.parse_args()
    backends = dict(item.split('=', 1) for item in arguments.backend)
    workers = [int(item) for item in arguments.workers.split(',')]
    if arguments.repetitions < 1 or any(value < 1 for value in workers):
        parser.error('Repetitions and worker counts must be positive')
    inputs = {str(path.resolve()): dict(bytes=path.stat().st_size, sha256=sha256(path)) for path in arguments.input}
    library_files = {label: {str(path.resolve()): sha256(path) for path in Path(directory).glob('*OpenMS*')
                            if path.is_file()} for label, directory in backends.items()}
    report = dict(platform=platform.platform(), machine=platform.machine(), python=sys.version,
                  cpu_count=os.cpu_count(), binary=str(arguments.binary.resolve()),
                  binary_sha256=sha256(arguments.binary), backends=backends, libraries=library_files, inputs=inputs,
                  parameters=vars(arguments) | {'binary': str(arguments.binary),
                    'input': [str(path) for path in arguments.input], 'output': str(arguments.output)},
                  records=[], summary=[])
    matrix = [(backend, filename, mode, worker, repetition)
              for backend in backends for filename in inputs
              for mode in arguments.modes.split(',') for worker in workers
              for repetition in range(arguments.repetitions)
              if mode == 'stream' or not filename.endswith(('.gz', '.bz2', '.zip'))]
    random.Random(20261019).shuffle(matrix)
    digests = {}
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    for number, (backend, filename, mode, worker, repetition) in enumerate(matrix, 1):
        environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        variable = 'DYLD_LIBRARY_PATH' if sys.platform == 'darwin' else 'LD_LIBRARY_PATH'
        environment[variable] = str(Path(backends[backend]).resolve()) + os.pathsep + environment.get(variable, '')
        command = [str(arguments.binary.resolve()), mode, filename, str(worker), str(arguments.queries)]
        row = measured(command, environment, arguments.timeout)
        if not row.get('openms_library') or Path(row['openms_library']).resolve().parent != Path(backends[backend]).resolve():
            raise RuntimeError(f'Wrong OpenMS library loaded for {backend}: {row.get("openms_library")}')
        row.update(backend=backend, input=filename, repetition=repetition)
        key = filename, mode
        fingerprint = (row['hash_per_worker'], row['spectra_per_worker'], row['chromatograms_per_worker'],
                       row['points_per_worker'], row['index_spectra'], row['index_chromatograms'])
        if key in digests and digests[key] != fingerprint:
            raise RuntimeError(f'Decoded-output mismatch for {key}: {digests[key]} != {fingerprint}, backend={backend}')
        digests[key] = fingerprint
        report['records'].append(row)
        report['summary'] = summarize(report['records'])
        arguments.output.write_text(json.dumps(report, indent=2) + '\n')
        print(f'{number}/{len(matrix)} {backend} {Path(filename).name} {mode} x{worker}: '
              f'{row["seconds"]:.3f}s, {row["max_rss_bytes"] / 1024 ** 2:.1f} MiB', flush=True)


if __name__ == '__main__':
    main()
