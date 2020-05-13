#! /bin/bash
# sh -> bash for debian. By J. R. Peterson. 2015/Jun.

pushd "`dirname "$0"`" > /dev/null 2>&1; rootdir="$PWD"; popd > /dev/null 2>&1;
MAFFT_BINARIES="$rootdir/mafftdir/libexec"; export MAFFT_BINARIES;

"$rootdir/mafftdir/bin/mafft" "$@"
# input file name can have space
