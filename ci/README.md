# Continuous integration

`.github/workflows/ci.yml` builds `proVBFH-inclusive`, `proVBFH`,
`proVBFHH` and `proVBFH-cs` on Linux and macOS against the latest releases
of hoppet, LHAPDF and FastJet, and runs short smoke tests, unit tests and
physics checks. It runs on every push and
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

`ci/nf-checks.sh` runs short non-factorisable (NF) runs of proVBFH and
proVBFHH with the gluon-mass regulator λ = M_V² and λ = 0.01 M_V²
(`nf_regfact`), which must give the same result: the log(λ/M_V²) terms
cancel point by point, and the check also tests the numerical angular
integrals against the analytic 1-loop triangle.

`proVBFH-cs` has no `configure`: it takes the tool paths and the sources
of its inclusive part from `../proVBFH`, so that is built first. Its CI
steps are the build, `make check` (unit tests of the kinematics, the
small-x channel and the four-parton generator), a smoke test of an NNLO
and an LO run, and the physics checks.

`ci/physics-checks.py` runs the physics checks (statistical comparisons
with pulls; see its `--help` for the options):

- `inclusive`: the inclusive cross section of `provbfh_incl`, `proVBFH`
  (`inclusive_only 1`) and the inclusive part of `proVBFH-cs` at LO, NLO,
  NNLO (and N3LO for the first two) and at the scale points (1,1), (1/2,1/2),
  (2,2), each against the frozen reference `ci/reference/inclusive.ref` and
  the implementations against each other. In `proVBFH-cs` the scale points
  are the on-the-fly ones of a single run, in the others separate runs.
- `projection`: in `proVBFH-cs` at NLO and NNLO the exclusive part must not
  change the Higgs-only observables (cross section without cuts, `ptH-incl`,
  `yH-incl`): they cancel point by point between the real radiation and its
  projection to Born, so these histograms vanish to rounding precision.
- `cs`: all histograms of `proVBFH-cs` at LO, NLO and NNLO and the three
  scale points against the frozen references `ci/reference/cs-<order>-W<k>.top`
  (chi2 over the bins, pull for the single-bin cross sections).

All of them use the hoppet benchmark PDF: `ci/benchmark-pdf.sh` evolves the
Les Houches benchmark initial condition with hoppet and writes it as a local
LHAPDF set (`hoppet_benchmark`, LHAPDF ID 900000, a few seconds), so the checks
do not depend on downloading a PDF set. The programs read it through the usual
LHAPDF route (`LHAPDF_DATA_PATH`).

The references are frozen: they come from long runs of the code at the
commit recorded in their headers, not from the run being tested, so a change
of the physics results (intended or not) is flagged. After an intended change,
regenerate them with `ci/physics-checks.py generate` (about an hour on 16
cores, mostly the NNLO exclusive part; `generate inclusive` or `generate
cs-nnlo` regenerate one part) and commit them with the change that needs them.
The sensitivity is limited by the statistics of the check runs: about 0.2%
on the inclusive cross sections, a few percent on the exclusive NLO and NNLO
cross sections with VBF cuts (more in the tails). `--stat 10` gives
sqrt(10) smaller errors at 10 times the run time, e.g. for a weekly run.

In the workflow each step runs through `ci/annotate.sh`, which reports
the error lines of a failing step as annotations in the run summary.

To check your current builds before pushing, run from the repository
root:

```
ci/run-local.sh                # all checks, on the programs as currently built
ci/run-local.sh proVBFHH       # only the checks of one package
ci/run-local.sh --build        # run make first
ci/run-local.sh --clean        # copy the working tree and build it from scratch, as CI does
ci/run-local.sh --help
```

By default nothing is copied or rebuilt: it warns if sources are newer
than an executable, and skips packages that are not built. Run
directories and logs go to a temporary directory. All checks take about
6 minutes.

The individual scripts can also be run directly, from the repository root:

```
ci/resolve-versions.sh                  # latest dependency versions
PREFIX=$HOME/ci-deps HOPPET_VERSION=... LHAPDF_VERSION=... FASTJET_VERSION=... \
    ci/install-deps.sh                  # build dependencies + PDF4LHC21_40
ci/smoke-tests.sh proVBFH-inclusive     # after building the package
ci/tensor-checks.sh                     # after building proVBFH-inclusive
ci/nf-checks.sh                         # after building proVBFH and proVBFHH
ci/physics-checks.py                    # after building proVBFH-inclusive, proVBFH, proVBFH-cs
ci/benchmark-pdf.sh DIR                 # the benchmark PDF alone, as an LHAPDF set in DIR
```
