#! /bin/sh

pushd "`dirname "$0"`" > /dev/null 2>&1; rootdir="$PWD"; popd > /dev/null 2>&1;
MAFFT_BINARIES="$rootdir/mafftdir/libexec"; export MAFFT_BINARIES;

"$rootdir/mafftdir/bin/mafft" "$@"
# $1 can have space in file name
