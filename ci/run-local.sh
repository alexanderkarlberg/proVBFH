#!/usr/bin/env bash
#
# Run the CI checks locally: smoke-tests.sh, tensor-checks.sh and
# nf-checks.sh, as in .github/workflows/ci.yml. By default they run on
# the programs as they are currently built in the repository, without
# copying or rebuilding anything. Continues after a failure and prints a
# summary; the exit status is non-zero if any check failed.
#
# Usage: ci/run-local.sh [options] [package ...]
#
#   package        proVBFH-inclusive, proVBFH and/or proVBFHH [default: all]
#                  (tensor-checks need proVBFH-inclusive, nf-checks
#                  proVBFH and/or proVBFHH)
#   (default)      use the existing builds in the repository; a package
#                  that is not built is skipped, and a warning is printed
#                  if sources are newer than its executable
#   --build        run make in each package first (and configure only if
#                  the package has never been configured, so an existing
#                  Makefile.inc is kept)
#   --clean        copy the working tree (uncommitted changes included)
#                  to the scratch directory and configure and build it
#                  from scratch, as CI does
#   --committed    the same for HEAD, i.e. what a push would test
#   --dir DIR      scratch directory for logs, run directories and copies
#                  [default: a new temporary directory]
#   --deps PREFIX  put PREFIX/bin first on PATH, e.g. an install made
#                  with ci/install-deps.sh with the CI's versions
#   -j N           parallel make jobs [default: number of cores]
#   -h, --help     this text
#
# The run directories always go to the scratch directory, never into the
# repository. The dependencies are found through their *-config programs
# on PATH, and the PDF set PDF4LHC21_40 must be installed.
#
set -uo pipefail

usage() { sed -n '/^# Usage/,/^# on PATH, and/p' "$0" | sed 's/^# \{0,1\}//'; }

MODE=existing
DIR=
JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
PKGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --build)     MODE=build ;;
        --clean)     MODE=clean ;;
        --committed) MODE=committed ;;
        --dir)       DIR=${2:?--dir needs an argument}; shift ;;
        --deps)      PATH="${2:?--deps needs an argument}/bin:$PATH"; shift ;;
        -j)          JOBS=${2:?-j needs an argument}; shift ;;
        -j*)         JOBS=${1#-j} ;;
        -h|--help)   usage; exit 0 ;;
        proVBFH-inclusive|proVBFH|proVBFHH) PKGS+=("$1") ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done
[ ${#PKGS[@]} -gt 0 ] || PKGS=(proVBFH-inclusive proVBFH proVBFHH)
has() { local p; for p in "${PKGS[@]}"; do [ "$p" = "$1" ] && return 0; done; return 1; }

REPO=$(git rev-parse --show-toplevel 2>/dev/null) || { echo "not in a git repository" >&2; exit 2; }
[ -n "$DIR" ] || DIR=$(mktemp -d "${TMPDIR:-/tmp}/provbfh-ci.XXXXXX")
mkdir -p "$DIR"
DIR=$(cd "$DIR" && pwd)
export RUNROOT="$DIR/runs"
LOGS="$DIR/logs"
mkdir -p "$RUNROOT" "$LOGS"

# executables of each package
exes() {
    case "$1" in
        proVBFH-inclusive) echo "proVBFH-inclusive/provbfh_incl proVBFH-inclusive/provbfhh_incl" ;;
        *) echo "$1/$1" ;;
    esac
}

# ----------------------------------------------------------------------
# the code to test
case "$MODE" in
    existing|build)
        SRC=$REPO ;;
    committed)
        SRC="$DIR/src"
        rm -rf "$SRC"; mkdir -p "$SRC"
        git -C "$REPO" archive HEAD | tar -x -C "$SRC" ;;
    clean)
        # tracked and new non-ignored files as they are on disk (deleted
        # files are skipped); releases/ and notes/ are not needed
        SRC="$DIR/src"
        rm -rf "$SRC"; mkdir -p "$SRC"
        (cd "$REPO" &&
         git ls-files -z --cached --others --exclude-standard -- . ':!releases' ':!notes' |
         while IFS= read -r -d '' f; do [ -e "$f" ] && printf '%s\0' "$f"; done |
         tar --null -T - -cf - | tar -x -C "$SRC") ;;
esac

changes=
[ "$MODE" != committed ] && [ -n "$(git -C "$REPO" status --porcelain --untracked-files=no -- . ':!releases' ':!notes')" ] &&
    changes=", with uncommitted changes"
echo "Mode:    $MODE ($(git -C "$REPO" rev-parse --abbrev-ref HEAD) at $(git -C "$REPO" rev-parse --short HEAD)$changes)"
echo "Code in: $SRC"
echo "Logs in: $LOGS"
for t in hoppet-config lhapdf-config fastjet-config; do
    if command -v $t > /dev/null; then
        echo "  $t $($t --version 2>/dev/null) ($(command -v $t))"
    else
        echo "  $t not found on PATH"
    fi
done
pdf_found=
IFS=: read -r -a pdfdirs <<< "${LHAPDF_DATA_PATH:-}:$(lhapdf-config --datadir 2>/dev/null)"
for d in "${pdfdirs[@]}"; do [ -n "$d" ] && [ -d "$d/PDF4LHC21_40" ] && pdf_found=1; done
[ -n "$pdf_found" ] || echo "WARNING: PDF set PDF4LHC21_40 not found; the runs will fail"
echo

# ----------------------------------------------------------------------
# steps
NAMES=(); RESULTS=()
FAILED=0

# step <name> <logfile> <command...>: run in $SRC, record the result
step() {
    local name=$1 log=$2 start=$SECONDS status
    shift 2
    printf '%-34s ' "$name"
    (cd "$SRC" && "$@") > "$log" 2>&1
    status=$?
    if [ $status -eq 0 ]; then
        RESULTS+=(ok); echo "ok      ($((SECONDS - start)) s)"
    else
        RESULTS+=(FAILED); FAILED=1
        echo "FAILED  ($((SECONDS - start)) s), last lines of $log:"
        tail -15 "$log" | sed 's/^/    | /'
    fi
    NAMES+=("$name")
    return $status
}
note() {
    printf '%-34s %s\n' "$1" "$2"
    NAMES+=("$1"); RESULTS+=("$2")
}

build() {
    local pkg=$1
    cd "$pkg" || return 1
    if [ "$MODE" = build ] && [ -f Makefile.inc ]; then
        make -j"$JOBS"
    else
        ./configure && make -j"$JOBS"
    fi
}

# is the package built? (prints why not; warns if its sources changed
# since the build)
built() {
    local pkg=$1 e newer
    for e in $(exes "$pkg"); do
        if [ ! -f "$SRC/$e" ]; then
            WHY="not built, skipped (build it, or use --build)"; return 1
        elif [ ! -x "$SRC/$e" ]; then
            # e.g. a sync (cernbox) that drops exec bits
            WHY="$e is not executable, skipped (chmod +x $e)"; return 1
        fi
        newer=$(cd "$SRC" && find "$pkg/src" "$pkg/analysis" "$pkg/Makefile" \
                    -newer "$e" -type f \( -name '*.f' -o -name '*.f90' -o -name '*.F' \
                    -o -name '*.h' -o -name '*.cc' -o -name Makefile \) 2>/dev/null | head -3)
        if [ -n "$newer" ]; then
            echo "WARNING: $e is older than" $newer "... (rebuild, or use --build)"
        fi
    done
    return 0
}

BUILT=" "   # packages available for the checks (bash 3 has no associative arrays)
for pkg in "${PKGS[@]}"; do
    if [ "$MODE" = existing ]; then
        if built "$pkg"; then
            BUILT="$BUILT$pkg "
        else
            note "$pkg" "$WHY"
            continue
        fi
    elif step "build $pkg" "$LOGS/build-$pkg.log" build "$pkg"; then
        BUILT="$BUILT$pkg "
    else
        continue
    fi
    step "smoke tests $pkg" "$LOGS/smoke-$pkg.log" ci/smoke-tests.sh "$pkg"
done

if [[ $BUILT == *" proVBFH-inclusive "* ]]; then
    # make check builds and runs the unit test (in the package's obj/)
    step "tensor checks" "$LOGS/tensor-checks.log" ci/tensor-checks.sh
fi

NF=()
for pkg in proVBFH proVBFHH; do [[ $BUILT == *" $pkg "* ]] && NF+=("$pkg"); done
if [ ${#NF[@]} -gt 0 ]; then
    step "NF checks (${NF[*]})" "$LOGS/nf-checks.log" ci/nf-checks.sh "${NF[@]}"
fi

# ----------------------------------------------------------------------
echo
if [ ${#RESULTS[@]} -eq 0 ] || ! printf '%s\n' "${RESULTS[@]}" | grep -qx ok; then
    echo "No checks ran."
    FAILED=1
elif [ $FAILED -eq 0 ]; then
    echo "All checks that ran passed."
    for i in "${!NAMES[@]}"; do
        [ "${RESULTS[$i]}" = ok ] || printf '  %-34s %s\n' "${NAMES[$i]}" "${RESULTS[$i]}"
    done
else
    echo "Some checks FAILED:"
    for i in "${!NAMES[@]}"; do
        [ "${RESULTS[$i]}" = ok ] || printf '  %-34s %s\n' "${NAMES[$i]}" "${RESULTS[$i]}"
    done
    echo "Logs: $LOGS, run directories: $RUNROOT"
fi
exit $FAILED
