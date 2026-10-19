#!/usr/bin/env python3
"""Verify complete Linux matrices and matching cross-mode decoded results."""
import json
from collections import defaultdict
from pathlib import Path
import sys
root=Path(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parents[1])
expected={'benchmarks.json':2560,'memory-scaling.json':960,'openmp-scaling.json':320,'metadata-scaling.json':640}
workers={1,2,4,8,16,32,64,96}
backends={'xerces','libxml2','expat','pugixml'}
reports={}
for name,count in expected.items():
 d=json.loads((root/name).read_text()); reports[name]=d
 assert not d['failures'],(name,d['failures'])
 assert len(d['records'])==count,(name,len(d['records']),count)
 groups=defaultdict(list)
 for r in d['records']:
  assert r['ok'],r
  w=r.get('openmp_threads',r['workers'])
  groups[r['backend'],Path(r['input']).name,r['mode'],w].append(r)
  assert r['max_rss_bytes']>0 and r['resident_rss_bytes']>0,r
 for key,rs in groups.items():
  assert len(rs)==5 and {r['repetition'] for r in rs}==set(range(5)),(name,key,rs)
  assert len({r['hash_per_worker'] for r in rs})==1,(name,key)
 assert {r['backend'] for r in d['records']}==backends,name
 assert {r.get('openmp_threads',r['workers']) for r in d['records']}==workers,name
 print(f'{name}: {count} records, {len(groups)} cells, five unique repetitions each')
stream={}; indexed={}; counts={}
for r in reports['benchmarks.json']['records']:
 name=Path(r['input']).name
 if r['mode']=='stream':
  key=name
  val=tuple(r[k] for k in ('hash_per_worker','spectra_per_worker','chromatograms_per_worker','points_per_worker'))
  assert key not in stream or stream[key]==val,key
  stream[key]=val
 elif r['mode']=='indexed':
  val=tuple(r[k] for k in ('hash_per_worker','spectra_per_worker','chromatograms_per_worker','points_per_worker'))
  assert name not in indexed or indexed[name]==val,name
  indexed[name]=val
 if r['mode'] in ('indexed','index-open'):
  val=(r['index_spectra'],r['index_chromatograms'])
  assert name not in counts or counts[name]==val,name
  counts[name]=val
for name,d in reports.items():
 for r in d['records']:
  file=Path(r['input']).name
  if r['mode']=='load' or name=='openmp-scaling.json':
   val=tuple(r[k] for k in ('hash_per_worker','spectra_per_worker','chromatograms_per_worker','points_per_worker'))
   assert val==stream[file],(name,r)
  elif r['mode']=='indexed-metadata':
   val=tuple(r[k] for k in ('hash_per_worker','spectra_per_worker','chromatograms_per_worker','points_per_worker'))
   assert val==indexed[file],(name,r)
  if r['mode'] in ('contexts','contexts-metadata','indexed-metadata'):
   assert (r['index_spectra'],r['index_chromatograms'])==counts[file],(name,r)
for name in ('benchmarks.json','memory-scaling.json','metadata-scaling.json'):
 assert reports[name]['libraries']==reports['benchmarks.json']['libraries'],name
print('All 4480 measurements complete; all counts, digests, cross-mode results and library hashes agree')
