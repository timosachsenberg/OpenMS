"""Fixed-work streaming: one reader, varying OpenMP binary-decoding threads.

Use a completed run_benchmarks.py report to select the identical binary,
backend snapshots, inputs and expected decoded results. Run without other jobs.
"""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import random
import statistics
import sys

from run_benchmarks import measured, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--threads', default='1,2,4,8')
    args = parser.parse_args()
    matrix = json.loads(args.matrix.read_text())
    thread_counts = [int(t) for t in args.threads.split(',')]
    assert all(t > 0 for t in thread_counts)
    binary = Path(matrix['binary'])
    launcher = binary.with_name('xml_benchmark_launcher')
    assert sha256(binary) == matrix['binary_sha256']
    assert sha256(launcher) == matrix['launcher_sha256']
    inputs = {p: data for p, data in matrix['inputs'].items()
              if Path(p).name in ('large.mzML', 'medium-zlib.mzML')}
    assert len(inputs) == 2
    assert all(sha256(p) == data['sha256'] for p, data in inputs.items())
    expected = {p: next(r['hash_per_worker'] for r in matrix['records']
                       if r['input'] == p and r['mode'] == 'stream') for p in inputs}
    report = {k:matrix[k] for k in ('platform','machine','binary','binary_sha256','launcher_sha256','backends','libraries')}
    report.update(kind='fixed-work OpenMP streaming', inputs=inputs, reader_workers=1,
                  openmp_threads=thread_counts, repetitions=5, records=[], failures=[], summary=[])
    jobs = [(b,p,t,r) for b in matrix['backends'] for p in inputs
            for t in thread_counts for r in range(5)]
    random.Random(20261020).shuffle(jobs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for n,(backend,filename,threads,repetition) in enumerate(jobs,1):
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1')
        variable = 'DYLD_LIBRARY_PATH' if sys.platform == 'darwin' else 'LD_LIBRARY_PATH'
        directory = Path(matrix['backends'][backend]).resolve()
        env[variable] = str(directory) + os.pathsep + env.get(variable,'')
        command = [str(launcher),str(binary),'stream',filename,'1','0']
        for attempt in range(2):
            row = measured(command,env,600)
            if row['ok']:
                break
            report['failures'].append(dict(row,backend=backend,input=filename,openmp_threads=threads,repetition=repetition,attempt=attempt))
            args.output.write_text(json.dumps(report,indent=2)+'\n')
            if row['phase'] != 'startup':
                raise RuntimeError(report['failures'][-1])
        if not row['ok']:
            raise RuntimeError(report['failures'][-1])
        assert Path(row['openms_library']).resolve().parent == directory
        assert row['hash_per_worker'] == expected[filename]
        row.update(backend=backend,input=filename,openmp_threads=threads,repetition=repetition)
        report['records'].append(row)
        groups=defaultdict(list)
        for item in report['records']:
            groups[item['backend'],item['input'],item['openmp_threads']].append(item)
        report['summary']=[dict(backend=b,input=p,openmp_threads=t,repetitions=len(rows),
            seconds_median=statistics.median(x['seconds'] for x in rows),
            seconds_min=min(x['seconds'] for x in rows),seconds_max=max(x['seconds'] for x in rows),
            rss_mib_median=statistics.median(x['max_rss_bytes']/1024**2 for x in rows))
            for (b,p,t),rows in sorted(groups.items())]
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(f'{n}/{len(jobs)} {backend} {Path(filename).name} OpenMP={threads}: {row["seconds"]:.3f}s',flush=True)


if __name__ == '__main__':
    main()
