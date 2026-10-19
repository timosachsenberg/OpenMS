#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
prefix="$root/dependencies/install"
mkdir -p "$root/source" "$root/staging" "$root/results" "$root/logs"
if ! test -f "$root/staging/openms-github.tar.gz"; then
 curl --fail --location --retry 3 https://codeload.github.com/OpenMS/OpenMS/tar.gz/77cccda7caf40b46a4919d33748450270eccd565 -o "$root/staging/openms-github.tar.gz"
fi
if test -d "$root/source/.git"; then echo "Use a fresh scratch directory; a source checkout already exists"; exit 1; fi
tar xzf "$root/staging/openms-github.tar.gz" --strip-components=1 -C "$root/source"
python3 "$root/staging/verify-source.py"
cd "$root/source"
if ! test -d .git; then
 git init -q
 printf '* -text\n' > .git/info/attributes
 git add -A
 git -c user.name='OpenMS XML POC' -c user.email=poc@localhost commit -qm 'Import exact archived OpenMS bytes'
 git tag imported-base
fi
mkdir -p "$root/dependencies/header-priority"
cd "$root/dependencies/header-priority"
for entry in xercesc expat.h expat_external.h pugixml.hpp pugiconfig.hpp zstd.h zstd_errors.h; do
 ln -sfn "$prefix/include/$entry" "$entry"
done
ln -sfn "$prefix/include/libxml2/libxml" libxml
printf '#include <xercesc/util/XercesVersion.hpp>\nstatic_assert(XERCES_VERSION_MAJOR == 3 && XERCES_VERSION_MINOR == 3);\n' | g++ -I"$root/dependencies/header-priority" -isystem /scratch/kohlbach/contrib/build/include -isystem "$prefix/include" -std=c++23 -fsyntax-only -x c++ -
