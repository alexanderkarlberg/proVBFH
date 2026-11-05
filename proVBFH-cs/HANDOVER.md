# proVBFH-cs: handover

State at commit `8008b27` on `2026-10-cs-production` (pushed; CI green on Ubuntu and macOS).

## What is done

- `qcd_order` selects the order in all codes; `cs_order` is a deprecated alias.
  `inclusive_only 1` runs only the inclusive part at any order.
- `make` in `proVBFH-cs` configures `../proVBFH` on demand; a fresh checkout builds.
- Link fix for macOS (`-lc++`).
- Benchmark PDF: `ci/benchmark-pdf.sh` writes hoppet's Les Houches benchmark as a local
  LHAPDF set (ID 900000), so the physics checks download nothing.
- `ci/physics-checks.py` (runs in CI, about 3 min):
  - `inclusive`: provbfh_incl, proVBFH and the cs inclusive part at LO, NLO, NNLO, N3LO
    and three scale points, against frozen references and each other.
  - `projection`: the Higgs-only observables are unchanged by the exclusive part (to 1e-14).
  - `cs`: every histogram at LO, NLO, NNLO against frozen references. Each histogram's
    integral (bins weighted by width) is a pull; the shape is a chi2 only where bins have
    relative error below 15%.
- Frozen references in `ci/reference/`, regenerated from commit `2188603` (headers say so).
- `make check` (proVBFH-cs) runs the unit tests and `tests/test_limits`, a pointwise test
  of the collinear cancellation of the (2,0)/(0,2) real: 1764 singular cases, worst residual
  3e-3 at y = 1e-8.
- `ci/limit-checks.sh` (weekly and on demand, about 4 min): the four-quark interference limits
  of `validation/v4_4qint/limits_4qint.f90`, with pass/fail criteria. Passes locally.

## Open items, in the order I would do them

1. **Check the next weekly run**, the first with `limit-checks.sh`. Also confirm the scheduled
   run is green.
2. **Shape coverage at NNLO is weak.** At 500k points per run, 37 to 49 histograms per run have
   no bin below 15% relative error, so their shapes are not compared; only integrals are.
   Options: more statistics for the NNLO reference (costly), or a coarser binning for the
   comparison.
3. **Exclusive-part tests not yet written:**
   - initial-final collinear limits (a final parton collinear to the incoming one);
   - soft limits;
   - the H+3j virtual against the VBFNLO reference points (the earlier bug showed up there).
   `tests/test_limits.f90` covers final-final collinear limits only.
4. **Sensitivity limits.** The physics checks would not catch a 10% change of the
   renormalisation scale or `cs_4qint 0`. The exclusive NLO/NNLO checks are only sensitive to
   gross errors at the statistics we run.
5. **Release (not started, deliberately).** Decided: `proVBFH-cs` stays its own directory and
   executable. Still to do: a `proVBFH-cs` entry in `CLAUDE.md`'s overview; `ChangeLog`,
   `AUTHORS`, `INSTALL`; tag; snapshot in `releases/`.
6. **Housekeeping.** Untracked in `1506.02660-results/`: the bug-fix symlink directory, the
   `plots-new/` PNGs and `__pycache__`. Untracked build products (`proVBFH/proVBFH`,
   `proVBFH/aux/combine_runs`, `proVBFH-inclusive/provbfh*`) are not all in `.gitignore`.

## How to run things

    # build (from proVBFH-cs; configures proVBFH on demand)
    make
    make check                       # unit tests + tests/test_limits, about 1.5 min
    # physics checks, all of them (about 3 min, needs the three executables built)
    ci/physics-checks.py
    ci/physics-checks.py cs --seed 9 # one part, another seed
    ci/limit-checks.sh               # interference limits, about 4 min
    # regenerate the references (about 1 h on 16 cores; do it from a clean commit)
    git worktree add --detach /tmp/refgen HEAD && ... ci/physics-checks.py generate

## Where things are

- Physics checks and frozen references: `ci/physics-checks.py`, `ci/reference/`.
- Benchmark PDF: `ci/benchmark-pdf.sh`, `ci/benchmark-pdf/make_benchmark_pdf.f90`.
- Pointwise tests: `tests/test_limits.f90` (and the test hook `nlo2_group_test` in
  `src/cs_nlo2.f90`, which only the test calls).
- Four-quark limits: `validation/v4_4qint/limits_4qint.f90`, reference `ref.f90`.
- CI: `.github/workflows/ci.yml`, `ci/README.md`.

## Things to keep in mind

- The collinear test found no bug in the dipoles; its first failures were in the test
  (wrong sign, stale maps, a criterion that also caught integrable terms). The fixes are in
  the commit message of `7804ff2`.
- The `cs_4qint` flag and the scale setting are documented in `proVBFH-cs/README`.
