# Continuous integration

`.github/workflows/ci.yml` builds `proVBFH-inclusive`, `proVBFH` and
`proVBFHH` on Linux and macOS against the latest releases of hoppet,
LHAPDF and FastJet, and runs short smoke tests. It runs on every push and
pull request, weekly (Mondays 05:00 UTC) to catch new dependency releases,
and on demand from the Actions tab, where specific dependency versions can
be pinned.

The smoke tests only check that the programs run to completion and write
finite results. The statistics are far too low for the numbers to mean
anything physically.

`ci/tensor-checks.sh` checks the tensor code: the unit test of
`tensor.f90` (`make check` in proVBFH-inclusive), that the three copies
of `tensor.f90` are identical, and that with zero widths `provbfhh_incl`
gives the same LO and NNLO cross section with the tensor (`-tensorME`)
and the analytic matrix element.

In the workflow each step runs through `ci/annotate.sh`, which reports
the error lines of a failing step as annotations in the run summary.

The scripts can also be run locally, from the repository root:

```
ci/resolve-versions.sh                  # latest dependency versions
PREFIX=$HOME/ci-deps HOPPET_VERSION=... LHAPDF_VERSION=... FASTJET_VERSION=... \
    ci/install-deps.sh                  # build dependencies + PDF4LHC21_40
ci/smoke-tests.sh proVBFH-inclusive     # after building the package
ci/tensor-checks.sh                     # after building proVBFH-inclusive
```
