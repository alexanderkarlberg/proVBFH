# Continuous integration

`.github/workflows/ci.yml` builds `proVBFH-inclusive`, `proVBFH`,
`proVBFHH` and `proVBFH-cs` on Linux and macOS against the latest releases
of hoppet, LHAPDF and FastJet, and runs short smoke tests, unit tests and
physics checks. It runs on every push and
pull request, weekly (Mondays 05:00 UTC) to catch new dependency releases,
and on demand from the Actions tab, where specific dependency versions can
be pinned.

`ci/limit-checks.sh` (weekly and on demand, about 4 minutes) runs the
single-unresolved limits of the identical-quark interference of the
four-quark reals of `proVBFH-cs` (`validation/v4_4qint`): the scalings of
the interference in the collinear and soft limits, the Cauchy-Schwarz bound,
and agreement of the code with N + I/2 in the limits.

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
with pulls, and exact comparisons of fixed runs; see its `--help` for the
options):

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
- `cs`: all histograms of `proVBFH-cs` at LO and the three
  scale points against the frozen references `ci/reference/cs-lo-W<k>.top`:
  the integral of each histogram (the bins weighted by their widths) is a pull,
  and the shape is a chi2 over the bins whose relative error is below 15%.
  At NLO and NNLO the `exact` check replaces this statistical comparison.
  (The references were generated with `ncall` 1e7, when `proVBFH-cs` mapped
  it to 1e6 for the inclusive part, the only part at LO; the checks and the
  generator now use 1e6 directly, which is the same statistics.)
- `exact` (Linux only, as in CI): `proVBFH-cs` at NLO and NNLO in short,
  fixed runs (fixed seed, fixed ncall) against the frozen outputs in
  `ci/reference/exact/`. The cross sections (`provbfh-*.dat`) must be reproduced
  byte for byte, and the histograms bin by bin (third column): the inclusive
  part (`provbfh-*-incl-*.top`, run with `cs_write_parts 1`) to 1e-5, the sum
  (`provbfh-*-W*.top`) and the exclusive part (`provbfh-*-excl-*.top`) to 1% of
  the bin or 3e-5 of the VBF cross section. The
  tolerance is not laxness: the runs are deterministic, but the exclusive
  part of `proVBFH-cs` contains discrete decisions (the radiation cutoff,
  the VBF cuts on the dipoles), so a last-bit difference of the phase space
  on another CPU moves a few events across a threshold and shifts a bin by
  the weight of those events. Measured across two machines with the same
  compiler and libraries, every inclusive observable is identical and the
  largest difference is 0.2% of a bin; the tolerance is about five times
  that, and still two orders of magnitude below any change of the physics.
  Any structural change (bins, histograms, number of files) still fails. A
  deliberate change of the output requires the regeneration
  `ci/physics-checks.py generate exact`, and the environment of the frozen
  outputs is recorded in `ci/reference/exact/README.md`. The frozen outputs
  are only comparable on the platform they were generated on, so the check
  runs on Linux only.

All of them use the hoppet benchmark PDF: `ci/benchmark-pdf.sh` evolves the
Les Houches benchmark initial condition with hoppet and writes it as a local
LHAPDF set (`hoppet_benchmark`, LHAPDF ID 900000, a few seconds), so the checks
do not depend on downloading a PDF set. The programs read it through the usual
LHAPDF route (`LHAPDF_DATA_PATH`).

The references are frozen: they come from long runs of the code at the
commit recorded in their headers, not from the run being tested, so a change
of the physics results (intended or not) is flagged. After an intended change,
regenerate them with `ci/physics-checks.py generate` (about an hour on 16
cores, mostly the NNLO exclusive part; `generate inclusive`, `generate cs-lo`
or `generate exact` regenerate one part) and commit them with the change that
needs them. The sensitivity of the statistical checks is limited by the
statistics of the check runs: about 0.2% on the inclusive cross sections.
`exact` is limited by the deterministic cross sections and the few-event
flips of the exclusive part, not by statistics (see above).

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
