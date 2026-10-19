"""Plot measured original/optimized medians; source-only input, no interpolation."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from plot_results import BACKENDS, LABELS, COLORS, cell


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',required=True,type=Path,help='Published results directory')
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for row,platform in enumerate(['macos','linux']):
        n=8 if platform=='macos' else 96
        groups={}
        for variant,directory in [('original',args.results/platform),('optimized',args.results/'optimized'/platform)]:
            groups[variant]={name:json.loads((directory/(name+'.json')).read_text())
                             for name in ['benchmarks','memory-scaling','metadata-scaling']}
            assert all(not d['failures'] and all(s['repetitions']==5 for s in d['summary']) for d in groups[variant].values())
        for col,(matrix,mode,title,unit,threads,metric,scale) in enumerate([
            ('benchmarks','stream','One-reader streaming','seconds; lower is better',1,'seconds_median',1),
            ('metadata-scaling','indexed-metadata','Default indexed open','seconds; lower is better',1,'job_open_seconds_median',1),
            ('memory-scaling','load',f'Whole-load peak RSS · {n} readers','GiB; lower is better',n,'rss_mib_median',1/1024)]):
            ax=axes[row,col]
            for b,backend in enumerate(BACKENDS):
                for j,variant in enumerate(['original','optimized']):
                    value=cell(groups[variant][matrix],backend,mode,threads)[metric]*scale
                    ax.bar(b+(-.19 if j==0 else .19),value,width=.36,color=COLORS[b],
                           alpha=.4 if j==0 else 1,hatch='//' if j==0 else None,
                           edgecolor=COLORS[b],linewidth=.6)
                    ax.annotate(f'{value:.2f}',(b+(-.19 if j==0 else .19),value),
                                xytext=(0,3 if j==0 else 12),textcoords='offset points',ha='center',fontsize=7)
            ax.set_xticks(range(4),LABELS);ax.set_ylabel(unit);ax.set_title(('macOS' if platform=='macos' else 'Linux')+' · '+title)
            ax.set_ylim(0,ax.get_ylim()[1]*1.15);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
            ax.spines[['right','top']].set_visible(False)
    fig.suptitle('OpenMS 3.7 mzML — original versus reviewed implementations',fontsize=17,y=.995)
    fig.legend(handles=[Patch(facecolor='#777',alpha=.4,hatch='//',label='Original frozen snapshot'),
                        Patch(facecolor='#777',label='Reviewed implementation')],
               loc='upper center',bbox_to_anchor=(.5,.96),ncol=2,frameon=False)
    fig.text(.5,.012,'Medians of 5 runs; identical ~1.05 GiB input, cached reads, decoding and digest included.\n'
             'Linux: 64 physical cores + 32 SMT siblings in the 96-reader mask. Shared fixes also applied to the Xerces control.\n'
             'Before/after sessions differ; small changes can reflect runtime or binary-layout variation. See raw ranges and control calibration.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.91))
    for suffix in ['png','svg']:fig.savefig(args.output/('optimization-comparison.'+suffix),dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':main()
