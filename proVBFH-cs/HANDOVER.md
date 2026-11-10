# proVBFH-cs: handover

State at commit `e85710b` on `2026-10-cs-production` (pushed; CI green on Ubuntu 26.04
and macOS, run 42828030165).

## What is done

- `qcd_order` selects the order in all codes; `cs_order` is a deprecated alias.
  `inclusive_only 1` runs only the inclusive part at any order.
- `make` in `proVBFH-cs` configures `../proVBFH` on demand; a fresh checkout builds.
- Link fix for macOS (`-lc++`).
- Benchmark PDF: `ci/benchmark-pdf.sh` writes hoppet's Les Houches benchmark as a local
  LHAPDF set (ID 900000), so the physics checks download nothing.
- `ci/physics-checks.py` (in CI about 1.5 min, plus 10 s for `exact` on Linux):
  - `inclusive`: provbfh_incl, proVBFH and the cs inclusive part at LO, NLO, NNLO, N3LO
    and three scale points, against frozen references and each other.
  - `projection`: the Higgs-only observables are unchanged by the exclusive part (to 1e-14).
  - `cs`: every histogram of proVBFH-cs at LO against frozen references. Each histogram's
    integral (bins weighted by width) is a pull; the shape is a chi2 only where bins have
    relative error below 15%. At NLO and NNLO `exact` replaces this statistical comparison.
  - `exact` (Linux only): proVBFH-cs at NLO and NNLO in short fixed runs (iseed 1,
    ncall1 = ncall2 = 2000, itmx1 = itmx2 = 1, 13 TeV, `cs_scales 3`, default card
    otherwise) against the frozen outputs in `ci/reference/exact/`. The cross sections
    (`xsct-*.dat`) must be reproduced byte for byte, the histograms bin by bin: the LO
    ones (`pwg-LO-*.top`) to 1e-5, the NLO/NNLO ones (`pwg-W*.top`, `pwg-EXCL-*.top`) to
    1% of the bin or 3e-5 of the VBF cross section. Histogram order, bin counts and the
    set of files are compared exactly. See "Things to keep in mind" for why the
    histograms need a tolerance.
- Frozen references: `ci/reference/inclusive.ref` and `ci/reference/cs-lo-W<k>.top` from
  commit `2188603` (2026-11-05), `ci/reference/exact/` from `9242a8b` (2026-11-10, the
  environment is recorded in its `README.md`).
- `make check` (proVBFH-cs) runs the unit tests and `tests/test_limits`, a pointwise test
  of the collinear cancellation of the (2,0)/(0,2) real: 1764 singular cases, worst residual
  3e-3 at y = 1e-8.
- `ci/limit-checks.sh` (weekly and on demand, about 4 min): the four-quark interference limits
  of `validation/v4_4qint/limits_4qint.f90`, with pass/fail criteria. Passes locally.
- CI: Linux is pinned to Ubuntu 26.04 with GCC 15.2 (`ubuntu-latest` is still 24.04), because
  the exact references are tied to the compiler; FastJet, hoppet and LHAPDF are built from
  the latest releases and cached. `ci/run-local.sh` runs the same steps locally (about 6 min).

## Open items, in the order I would do them

1. **Check the next weekly run**, the first with `limit-checks.sh` in the schedule.
   Also confirm the scheduled run is green.
2. **The exact check pins the NLO/NNLO histograms to 1%, not digit for digit.** Bins that
   are empty in the 2000-point reference are only bounded to 3e-5 of the VBF cross section
   (about 3e-5 pb). Tighter would need more points in the exact runs (still seconds) or a
   smaller cross-machine spread; measure the spread before tightening, not after.
3. **Exclusive-part tests not yet written:**
   - initial-final collinear limits (a final parton collinear to the incoming one);
   - soft limits;
   - the H+3j virtual against the VBFNLO reference points (the earlier bug showed up there).
   `tests/test_limits.f90` covers final-final collinear limits only.
4. **Non-default options are not covered.** The exact runs use the whole default card
   (`cs_scales 3`, `cs_smallx 0.3`, `cs_4qint 1`), so any change of the integrand there moves
   `xsct` and fails. Not exercised: other scale choices, `cs_4qint 0`, other `cs_smallx`
   values, `inclusive_only`, and the `cs_4qint 2` subtraction (only in the validation).
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
    # physics checks, all of them (about 1 min; about 1.5 min in CI,
    # needs the three executables built)
    ci/physics-checks.py
    ci/physics-checks.py exact       # only the NLO/NNLO fixed runs, about 10 s
    ci/physics-checks.py cs --seed 9 # one part, another seed
    ci/limit-checks.sh               # interference limits, about 4 min
    ci/run-local.sh                  # everything CI runs, about 6 min
    # regenerate the references (about 1 h on 16 cores; do it from a clean commit)
    git worktree add --detach /tmp/refgen HEAD && ... ci/physics-checks.py generate
    ci/physics-checks.py generate exact   # only the exact references, about 10 s

## Where things are

- Physics checks and frozen references: `ci/physics-checks.py`, `ci/reference/`.
- Benchmark PDF: `ci/benchmark-pdf.sh`, `ci/benchmark-pdf/make_benchmark_pdf.f90`.
- Pointwise tests: `tests/test_limits.f90` (and the test hook `nlo2_group_test` in
  `src/cs_nlo2.f90`, which only the test calls).
- Four-quark limits: `validation/v4_4qint/limits_4qint.f90`, reference `ref.f90`.
- CI: `.github/workflows/ci.yml`, `ci/README.md`, `ci/run-local.sh`.

## Things to keep in mind

- Byte-identical *histograms* across machines are not attainable and not worth chasing:
  the exclusive part makes discrete decisions (the radiation cutoff, the VBF cuts on the
  dipoles), so a last-bit phase-space difference moves a few events across a threshold.
  Between an Intel desktop and the runner (same compiler, hoppet and LHAPDF; FastJet 3.5.3
  was also checked locally and changes nothing) every inclusive observable and `xsct` are
  identical, and the worst bin differs by 0.2% of the bin (1.2e-4 of the VBF cross
  section); the tolerance is about five times that. Ruled out as causes: OpenMP,
  FastJet 3.5.1 vs 3.5.3, the PDF grid, uninitialised heap, and
  `GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-AVX512F` (tried, removed: it moves the
  results by 1e-8 but does not bridge the machines).
- Ubuntu 26.04 ships GCC 15.2.0, which is why the runner is pinned to it: the exact
  references come from that compiler. A compiler or library update means regenerating
  `ci/reference/exact/` deliberately, not loosening the check.
- The collinear test found no bug in the dipoles; its first failures were in the test
  (wrong sign, stale maps, a criterion that also caught integrable terms). The fixes are in
  the commit message of `7804ff2`.
- The `cs_4qint` flag and the scale setting are documented in `proVBFH-cs/README`.
