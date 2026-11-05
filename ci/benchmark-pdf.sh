#!/usr/bin/env bash
#
# Generates the hoppet benchmark PDF (Les Houches benchmark initial
# condition evolved with hoppet) as a local LHAPDF set, so that the
# physics checks do not need to download a PDF set from the LHAPDF server.
#
# Usage: ci/benchmark-pdf.sh <dir>
#
# creates <dir>/hoppet_benchmark/{hoppet_benchmark.info,hoppet_benchmark_0000.dat}
# and <dir>/pdfsets.index (LHAPDF ID 900000), and prints <dir>. Add the
# directory to the LHAPDF search path with
#     export LHAPDF_DATA_PATH=<dir>:$(lhapdf-config --datadir)
# (the checks do this through ci/common.sh). Needs hoppet (>= 2.1) and
# gfortran; takes a few seconds. An existing set in <dir> is reused.
#
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
DIR=${1:?usage: $0 <dir>}
NAME=hoppet_benchmark
ID=900000

mkdir -p "$DIR"
DIR=$(cd "$DIR" && pwd)
if [ ! -s "$DIR/$NAME/${NAME}_0000.dat" ] || [ ! -s "$DIR/pdfsets.index" ]; then
    work=$(mktemp -d)
    trap 'rm -rf "$work"' EXIT
    ${FC:-gfortran} $(hoppet-config --fflags) -J"$work" "$HERE/benchmark-pdf/make_benchmark_pdf.f90" \
        $(hoppet-config --libs) -o "$work/make_benchmark_pdf"
    rm -rf "$DIR/$NAME"
    mkdir -p "$DIR/$NAME"
    "$work/make_benchmark_pdf" "$DIR/$NAME/$NAME" > "$work/log" 2>&1 || { cat "$work/log" >&2; exit 1; }
    sed -i.bak -e "s/^SetIndex:.*/SetIndex: $ID/" \
               -e "s/^SetDesc:.*/SetDesc: \"Les Houches benchmark PDF evolved with hoppet\"/" \
               -e "s/^Authors:.*/Authors: hoppet/" \
               -e "s/^Reference:.*/Reference: \"hoppet-code\/hoppet\"/" \
               -e "s/^ErrorType:.*/ErrorType: null/" "$DIR/$NAME/$NAME.info"
    rm -f "$DIR/$NAME/$NAME.info.bak"
    echo "$ID $NAME 1" > "$DIR/pdfsets.index"
fi
echo "$DIR"
