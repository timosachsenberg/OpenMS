#!/bin/bash
set -euo pipefail
original=${XML_POC_ROOT}
root="$original/optimized"
mkdir -p "$root"/{source,staging,logs,results,snapshots,tmp}
ln -s "$original/dependencies" "$root/dependencies"
ln -s "$original/data" "$root/data"
cd "$root/source"
git init -q
printf '* -text\n' > .git/info/attributes
git fetch --depth 1 https://github.com/OpenMS/OpenMS.git 77cccda7caf40b46a4919d33748450270eccd565
git fetch "$root/staging/source.bundle" '+refs/heads/*:refs/heads/*'
git bundle verify "$root/staging/source.bundle"
