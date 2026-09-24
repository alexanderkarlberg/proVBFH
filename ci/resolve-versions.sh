#!/usr/bin/env bash
#
# Print the latest released versions of the external dependencies, one
# KEY=VALUE per line (suitable for appending to $GITHUB_OUTPUT).
#
# Any version can be pinned by setting the corresponding environment
# variable (HOPPET_VERSION, LHAPDF_VERSION, FASTJET_VERSION) beforehand.
#
set -euo pipefail

# GitHub API (authenticated if GH_TOKEN is set, to avoid rate limits)
auth=()
if [ -n "${GH_TOKEN:-}" ]; then auth=(-H "Authorization: Bearer $GH_TOKEN"); fi

if [ -z "${HOPPET_VERSION:-}" ]; then
    HOPPET_VERSION=$(curl -fsSL "${auth[@]}" \
        https://api.github.com/repos/hoppet-code/hoppet/releases/latest \
        | sed -n 's/.*"tag_name": *"hoppet-\([0-9.]*\)".*/\1/p')
fi
if [ -z "${LHAPDF_VERSION:-}" ]; then
    LHAPDF_VERSION=$(curl -fsSL https://lhapdf.hepforge.org/downloads \
        | grep -oE 'LHAPDF-6\.[0-9]+\.[0-9]+\.tar\.gz' \
        | sed 's/LHAPDF-//; s/\.tar\.gz//' | sort -uV | tail -1)
fi
if [ -z "${FASTJET_VERSION:-}" ]; then
    FASTJET_VERSION=$(curl -fsSL https://fastjet.fr/all-releases.html \
        | grep -oE 'fastjet-3\.[0-9]+\.[0-9]+\.tar\.gz' \
        | sed 's/fastjet-//; s/\.tar\.gz//' | sort -uV | tail -1)
fi

for v in HOPPET_VERSION LHAPDF_VERSION FASTJET_VERSION; do
    if [ -z "${!v}" ]; then
        echo "ERROR: could not determine $v" >&2
        exit 1
    fi
done

echo "hoppet=$HOPPET_VERSION"
echo "lhapdf=$LHAPDF_VERSION"
echo "fastjet=$FASTJET_VERSION"
