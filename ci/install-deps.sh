#!/usr/bin/env bash
#
# Build LHAPDF, FastJet and hoppet from their release tarballs and
# install them into $PREFIX, then download the PDF sets in $PDFSETS.
#
# Required environment:
#   PREFIX                                        install location
#   HOPPET_VERSION, LHAPDF_VERSION, FASTJET_VERSION  (see resolve-versions.sh)
# Optional:
#   PDFSETS   space-separated LHAPDF set names [default: PDF4LHC21_40]
#   FC        Fortran compiler [default: gfortran]
#
set -euo pipefail

: "${PREFIX:?}" "${HOPPET_VERSION:?}" "${LHAPDF_VERSION:?}" "${FASTJET_VERSION:?}"
PDFSETS=${PDFSETS:-PDF4LHC21_40}
FC=${FC:-gfortran}
NCPU=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)

mkdir -p "$PREFIX"
PREFIX=$(cd "$PREFIX" && pwd)
export PATH="$PREFIX/bin:$PATH"
BUILD=$(mktemp -d)
cd "$BUILD"

echo "::group::LHAPDF $LHAPDF_VERSION"
curl -fsSL "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz" | tar xz
(cd "LHAPDF-$LHAPDF_VERSION" \
     && ./configure --prefix="$PREFIX" --disable-python \
     && make -j"$NCPU" && make install)
echo "::endgroup::"

echo "::group::FastJet $FASTJET_VERSION"
curl -fsSL "https://fastjet.fr/repo/fastjet-$FASTJET_VERSION.tar.gz" | tar xz
(cd "fastjet-$FASTJET_VERSION" \
     && ./configure --prefix="$PREFIX" \
     && make -j"$NCPU" && make install)
echo "::endgroup::"

echo "::group::hoppet $HOPPET_VERSION"
curl -fsSL "https://github.com/hoppet-code/hoppet/archive/refs/tags/hoppet-$HOPPET_VERSION.tar.gz" | tar xz
cmake -S "hoppet-hoppet-$HOPPET_VERSION" -B hoppet-build \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_Fortran_COMPILER="$FC" \
      -DHOPPET_BUILD_EXAMPLES=OFF \
      -DHOPPET_ENABLE_TESTING=OFF \
      -DHOPPET_BUILD_BENCHMARK=OFF \
      -DHOPPET_GIT_WATCHER=OFF
cmake --build hoppet-build -j "$NCPU"
cmake --install hoppet-build
echo "::endgroup::"

echo "::group::PDF sets: $PDFSETS"
DATADIR=$(lhapdf-config --datadir)
for set in $PDFSETS; do
    curl -fsSL "https://lhapdfsets.web.cern.ch/current/$set.tar.gz" | tar xz -C "$DATADIR"
done
echo "::endgroup::"

cd /
rm -rf "$BUILD"

hoppet-config --version
lhapdf-config --version
fastjet-config --version
