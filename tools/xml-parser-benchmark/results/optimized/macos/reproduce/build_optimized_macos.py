"""Build isolated optimized snapshots without changing the frozen original experiments."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'OpenMS-xerces-control'
BUILD=ROOT/'OpenMS-optimized-build'
BENCH=ROOT/'benchmark-optimized-build'
OUTPUT=ROOT/'results/optimized/macos'
OUTPUT.mkdir(parents=True,exist_ok=True)
EXISTING=('OnDiscMSExperiment ConsensusXMLFile FeatureXMLFile IndexedMzMLDecoder IndexedMzMLFile '
          'IndexedMzMLFileLoader MzMLFile MzMLSpectrumDecoder TransformationXMLFile XMLFile XMLHandler '
          'XMLValidator CVMappingFile IdXMLFile MzDataFile MzIdentMLFile MzDataValidator MzMLValidator '
          'MzXMLFile ParamXMLFile PepXMLFile ProtXMLFile SemanticValidator TraMLFile UnimodXMLFile').split()
FOCUSED={'libxml2':'XMLParserBackend_test','expat':'ExpatXMLFile_test','pugixml':'XMLBackendPrototype_test'}
REFS={'xerces':'experiment/xml-xerces-control','libxml2':'feature/xml-libxml2-poc',
      'expat':'feature/xml-expat-poc','pugixml':'feature/xml-pugixml-poc'}


def run(args,log,env=None):
    print(' '.join(str(a) for a in args[:5]),flush=True)
    with (OUTPUT/log).open('w') as handle:
        subprocess.run(list(map(str,args)),cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT,check=True,env=env)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--backends',default=','.join(REFS))
    args=parser.parse_args()
    manifest={}
    path=OUTPUT/'source-manifest.json'
    if path.exists():manifest=json.loads(path.read_text())
    for backend in args.backends.split(','):
        revision=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse',REFS[backend]],text=True).strip()
        run(['git','-C',SOURCE,'checkout','--detach',revision],backend+'-source.log')
        run(['cmake','-S',SOURCE,'-B',BUILD,'-G','Ninja','-DCMAKE_BUILD_TYPE=Release',
             '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
             '-DOPENMS_USE_VCPKG=OFF','-DWITH_GUI=OFF','-DBUILD_TOPP_TOOLS=OFF','-DPYOPENMS=OFF',
             '-DENABLE_DOCS=OFF','-DENABLE_CLASS_TESTING=ON','-DENABLE_TOPP_TESTING=OFF',
             '-DWITH_OPENTIMS=OFF','-DWITH_THERMO_RAW=OFF','-DWITH_ONNX=OFF','-DWITH_WNETALIGN=OFF',
             '-DGIT_TRACKING=OFF','-DENABLE_UPDATE_CHECK=OFF','-DMT_REQUIRE_OPENMP=ON',
             '-DCMAKE_PREFIX_PATH=/opt/homebrew;/opt/homebrew/opt/libxml2;/opt/homebrew/opt/libomp',
             '-DCURL_INCLUDE_DIR=/opt/homebrew/opt/curl/include',
             '-DCURL_LIBRARY=/opt/homebrew/opt/curl/lib/libcurl.dylib',
             '-DCURL_LIBRARY_RELEASE=/opt/homebrew/opt/curl/lib/libcurl.dylib',
             '-DLIBXML2_INCLUDE_DIR=/opt/homebrew/opt/libxml2/include/libxml2',
             '-DLIBXML2_LIBRARY=/opt/homebrew/opt/libxml2/lib/libxml2.dylib',
             '-DEXPAT_INCLUDE_DIR='+str(ROOT/'dependencies/expat-install/include'),
             '-DEXPAT_LIBRARY='+str(ROOT/'dependencies/expat-install/lib/libexpat.dylib')],backend+'-configure.log')
        targets=[name+'_test' for name in EXISTING]
        if backend in FOCUSED:targets.append(FOCUSED[backend])
        run(['cmake','--build',BUILD,'--target',*targets,'-j8'],backend+'-build.log')
        run(['ctest','--test-dir',BUILD,'-R','^('+'|'.join(targets)+')$',
             '--output-on-failure','-j1'],backend+'-tests.log')
        snapshot=ROOT/'snapshots-optimized'/backend/'lib'
        snapshot.mkdir(parents=True,exist_ok=True)
        for pattern in ('libOpenMS*.dylib','libOpenSwathAlgo*.dylib'):
            for library in (BUILD/'lib').glob(pattern):shutil.copy2(library,snapshot/library.name)
        if backend=='expat':
            for library in (ROOT/'dependencies/expat-install/lib').glob('libexpat*.dylib'):
                shutil.copy2(library,snapshot/library.name)
        shutil.copy2(BUILD/'compile_commands.json',OUTPUT/(backend+'-compile_commands.json'))
        shutil.copy2(BUILD/'CMakeCache.txt',OUTPUT/(backend+'-build-cache.txt'))
        source_tree=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse',revision+'^{tree}'],text=True).strip()
        manifest[backend]={'source_commit':revision,'source_tree':source_tree,'optimization':'Release -O3 -DNDEBUG, C++23',
                           'timed_source_sha256':{name:hashlib.sha256((SOURCE/'tools/xml-parser-benchmark'/name).read_bytes()).hexdigest()
                              for name in ('benchmark.cpp','launcher.cpp','run_benchmarks.py','run_openmp_scaling.py')},
                           'libraries':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in snapshot.glob('*.dylib')}}
        path.write_text(json.dumps(manifest,indent=2)+'\n')
        if not (BENCH/'xml_parser_benchmark').exists():
            run(['cmake','-S',SOURCE/'tools/xml-parser-benchmark','-B',BENCH,'-G','Ninja',
                 '-DCMAKE_BUILD_TYPE=Release','-DOpenMS_DIR='+str(BUILD),
                 '-DOpenMP_CXX_FLAGS=-Xclang -fopenmp','-DOpenMP_CXX_LIB_NAMES=omp',
                 '-DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib',
                 '-DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include'], 'benchmark-configure.log')
            run(['cmake','--build',BENCH,'-j8'],'benchmark-build.log')
        import os
        environment=dict(os.environ,DYLD_LIBRARY_PATH=str(snapshot),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        run(['ctest','--test-dir',BENCH,'--output-on-failure'],backend+'-benchmark-tests.log',environment)
        verification=subprocess.run(list(map(str,[BENCH/'xml_benchmark_launcher',BENCH/'xml_parser_benchmark','verify',
             ROOT/'benchmark-data/huge-array.mzML','1','2048'])),cwd=ROOT,env=environment,capture_output=True,text=True,check=True)
        (OUTPUT/(backend+'-huge-array-stderr.log')).write_text(verification.stderr)
        result=json.loads(next(line for line in verification.stdout.splitlines() if line.startswith('{')))
        (OUTPUT/(backend+'-huge-array-equivalence.json')).write_text(json.dumps(result,indent=2)+'\n')
        print('COMPLETE',backend,revision,flush=True)


if __name__=='__main__':main()
