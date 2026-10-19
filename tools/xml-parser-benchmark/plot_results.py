"""Plot paired platform results; requires matplotlib (no data-frame dependency)."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

BACKENDS = ('xerces', 'libxml2', 'expat', 'pugixml')
LABELS = ('Xerces-C++', 'libxml2', 'Expat', 'pugixml')
COLORS = ('#4b5563', '#0072b2', '#009e73', '#d55e00')


def cell(data, backend, mode, workers, filename='large.mzML'):
    return next(r for r in data['summary'] if r['backend'] == backend and
                r.get('mode') == mode and r['workers'] == workers and
                Path(r['input']).name == filename)


def axes_style(ax, counts, linux=False):
    ax.set_xscale('log', base=2)
    ax.set_xticks(counts)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlim(.9, max(counts)*1.1)
    ax.set_xlabel('Threads / independent readers')
    ax.grid(alpha=.2, which='major')
    ax.spines[['right','top']].set_visible(False)
    if linux and max(counts) > 64:
        ax.axvspan(64, max(counts), color='grey', alpha=.10)
        ax.axvline(64, color='grey', linewidth=.8, linestyle=':')


def save(fig, name, title, footer, out):
    fig.suptitle(title, fontsize=16, y=.995)
    handles = [Line2D([], [], color=c, marker='o', label=l) for c,l in zip(COLORS,LABELS)]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5,.955), ncol=4, frameon=False)
    fig.text(.5,.012,footer,ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.055,1,.90))
    for suffix in ('png','svg'):
        fig.savefig(out/(name+'.'+suffix), dpi=170, bbox_inches='tight')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mac',required=True,type=Path)
    parser.add_argument('--linux',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    platforms=[]
    for label,directory in (('macOS · Apple M4 Max',args.mac),('Linux · AMD EPYC 7763',args.linux)):
        reports={k:json.loads((directory/(k+'.json')).read_text())
                 for k in ('benchmarks','memory-scaling','openmp-scaling','metadata-scaling')}
        assert all(not v['failures'] for v in reports.values())
        assert all(r['repetitions']==5 for v in reports.values() for r in v['summary'])
        residency = directory/'context-residency.json'
        reports['context-residency'] = json.loads(residency.read_text()) if residency.exists() else reports['memory-scaling']
        counts=sorted({r['workers'] for r in reports['benchmarks']['summary']})
        platforms.append((label,reports,counts))
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'svg.fonttype':'none'})
    modes=('stream','load','indexed-metadata')
    titles=('Streaming API / digest consumer','Retained whole-dataset load','Default indexed open + reads')
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for row,(label,reports,counts) in enumerate(platforms):
        for col,(mode,title) in enumerate(zip(modes,titles)):
            ax=axes[row,col]
            data=reports[{'load':'memory-scaling','indexed-metadata':'metadata-scaling'}.get(mode,'benchmarks')]
            ax.plot(counts,counts,':',color='#999999',linewidth=1,label='Ideal')
            for b,color in zip(BACKENDS,COLORS):
                one=cell(data,b,mode,1)['seconds_median']
                rows=[cell(data,b,mode,n) for n in counts]
                gains=[n*one/r['seconds_median'] for n,r in zip(counts,rows)]
                ax.plot(counts,gains,'o-',color=color,linewidth=1.6,markersize=4)
                ax.fill_between(counts,[n*one/r['seconds_max'] for n,r in zip(counts,rows)],
                                [n*one/r['seconds_min'] for n,r in zip(counts,rows)],color=color,alpha=.10)
            axes_style(ax,counts,row==1)
            ax.set_yscale('log',base=2)
            ax.set_ylim(.7,max(counts)*1.4)
            ax.set_yticks([n for n in (1,2,4,8,16,32,64,128) if n<=max(counts)*1.4])
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.set_ylabel('Aggregate throughput gain (×)')
            ax.set_title(label+'\n'+title)
    save(fig,'concurrent-scaling','Concurrent reader throughput — identical ~1 GiB mzML',
         'Median of 5; min/max bands clipped at plot limits (full ranges in CSV). Each curve is normalized to its own one-reader baseline. Each reader repeats the whole workload.\n'
         'Dotted line: ideal throughput. Linux shading: beyond 64 physical cores (SMT). Different hosts/toolchains; cached inputs.',args.output)
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for row,(label,reports,counts) in enumerate(platforms):
        for col,(mode,title) in enumerate(zip(modes,titles)):
            ax=axes[row,col]
            data=reports[{'load':'memory-scaling','indexed-metadata':'metadata-scaling'}.get(mode,'benchmarks')]
            for b,color in zip(BACKENDS,COLORS):
                ax.plot(counts,[cell(data,b,mode,n)['rss_mib_median']/1024 for n in counts],
                        'o-',color=color,linewidth=1.6,markersize=4)
            axes_style(ax,counts,row==1)
            ax.set_yscale('log',base=10)
            ax.set_ylabel('Whole-process peak RSS (GiB)')
            ax.set_title(label+'\n'+title)
    save(fig,'memory-scaling','Peak resident memory — streaming, retained loads, parallel random reads',
         'Median of 5 process peaks; includes runtime/CV initialization. Same ~1 GiB file. Pugixml retains XML DOM while parsing.\n'
         'Loaded experiments remain alive until all workers finish. Linux shading marks SMT beyond 64 physical cores.',args.output)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for row,(label,reports,counts) in enumerate(platforms):
        idle=reports['context-residency']
        for col,(report_name,mode,title) in enumerate((
                ('context-residency','contexts','Index only; metadata skipped'),
                ('metadata-scaling','contexts-metadata','Default open; metadata retained'))):
            ax=axes[row,col]; data=reports[report_name]
            for b,color in zip(BACKENDS,COLORS):
                values=[(cell(data,b,mode,n)['resident_rss_mib_median']-
                         cell(idle,b,'idle',n)['resident_rss_mib_median'])/n for n in counts]
                ax.plot(counts,[v if v>0 else float('nan') for v in values],'o-',color=color,markersize=4)
            axes_style(ax,counts,row==1)
            ax.set_yscale('log')
            ax.set_ylabel('Incremental resident RSS (MiB / reader)')
            ax.set_title(label+'\n'+title)
    save(fig,'reader-context-memory','Open-reader resident overhead — same ~1 GiB mzML',
         '(Median resident RSS with contexts − median resident RSS with idle workers) / readers. All opened contexts remain alive.\n'
         'Sampled after joining; includes indexes/metadata and allocator retention, not sizeof a native parser; Linux resident counters are approximate. Nonpositive noisy estimates omitted.',args.output)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for row,(label,reports,_) in enumerate(platforms):
        data=reports['openmp-scaling']; counts=data['openmp_threads']
        for col,filename in enumerate(('large.mzML','medium-zlib.mzML')):
            ax=axes[row,col]
            for b,color in zip(BACKENDS,COLORS):
                rows=sorted([r for r in data['summary'] if r['backend']==b and Path(r['input']).name==filename],key=lambda r:r['openmp_threads'])
                assert [r['openmp_threads'] for r in rows]==counts
                one=rows[0]['seconds_median']
                ax.plot(counts,[one/r['seconds_median'] for r in rows],'o-',color=color,markersize=4)
            axes_style(ax,counts,row==1)
            ax.set_xlabel('Requested OpenMP team; one reader, fixed input')
            ax.axhline(1,color='grey',linestyle=':',linewidth=1)
            ax.set_ylabel('Single-file speedup (×)')
            ax.set_title(label+'\n'+('~1 GiB plain mzML' if col==0 else '~66 MiB binary-zlib mzML'))
    save(fig,'openmp-scaling','Single-file binary-decoding parallelism — fixed work',
         'Median of 5, normalized to the same backend at OpenMP=1. Parsing and digest consumption remain serial; OpenMS binary-array decoding uses OpenMP.\n'
         'Includes decoded consumption. Linux shading marks SMT beyond 64 physical cores; runtime can retain extra initialization threads.',args.output)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for row,(label,reports,_) in enumerate(platforms):
        data=reports['openmp-scaling']; counts=data['openmp_threads']
        for col,filename in enumerate(('large.mzML','medium-zlib.mzML')):
            ax=axes[row,col]
            for b,color in zip(BACKENDS,COLORS):
                rows=sorted([r for r in data['summary'] if r['backend']==b and Path(r['input']).name==filename],key=lambda r:r['openmp_threads'])
                ax.plot(counts,[r['rss_mib_median']/1024 for r in rows],'o-',color=color,markersize=4)
            axes_style(ax,counts,row==1)
            ax.set_xlabel('Requested OpenMP team; one reader, fixed input')
            ax.set_yscale('log')
            ax.set_ylabel('Whole-process peak RSS (GiB)')
            ax.set_title(label+'\n'+('~1 GiB plain mzML' if col==0 else '~66 MiB binary-zlib mzML'))
    save(fig,'openmp-memory','Resident memory with parallel binary decoding — one reader',
         'Median of 5 process peaks. Same fixed input at each thread count; whole-process RSS includes initialization.\n'
         'Linux shading marks SMT beyond 64 physical cores. RSS includes retained initialization threads; process thread count can exceed the requested team size.',args.output)
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for col,(label,reports,counts) in enumerate(platforms):
        for row in (0,1):
            ax=axes[row,col]
            for b,color in zip(BACKENDS,COLORS):
                for report_name,mode,style in (('benchmarks','indexed','--'),('metadata-scaling','indexed-metadata','-')):
                    rows=[cell(reports[report_name],b,mode,n) for n in counts]
                    values=[2048*n/r['seconds_median'] for n,r in zip(counts,rows)] if row==0 else [r['rss_mib_median']/1024 for r in rows]
                    ax.plot(counts,values,color=color,linestyle=style,marker='o',markersize=3)
            axes_style(ax,counts,col==1)
            ax.set_yscale('log')
            ax.set_ylabel('Aggregate spectrum requests / second' if row==0 else 'Whole-process peak RSS (GiB)')
            ax.set_title(label)
    save(fig,'indexed-path-comparison','Indexed access — default metadata versus fast index-only open',
         'Solid: default metadata loading. Dashed: skipMetaData=true. 2,048 identical random requests per reader, same ~1 GiB file.\n'
         'Median of 5; total timing and peak RSS include opening. Linux shading marks SMT beyond 64 physical cores.',args.output)


if __name__=='__main__':
    main()
