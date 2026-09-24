# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository overview

This repository contains Fortran/POWHEG-BOX physics codes computing QCD corrections to Vector Boson
Fusion (VBF) Higgs and Higgs-pair production, developed for phenomenology papers (see `notes/`). It
contains three independent, self-contained programs, each with its own `configure`/`Makefile`,
`src/`, `docs/`, `AUTHORS`, `ChangeLog`, `README`, `INSTALL`:

- **`proVBFH/`** — fully differential single-Higgs VBF production, built on the POWHEG BOX framework
  (NLO+PS matching), extended to NNLO via the structure-function approach and to N3LO in inclusive mode.
- **`proVBFHH/`** — the double-Higgs (VBF HH) analogue of `proVBFH`, same architecture.
- **`proVBFH-inclusive/`** — a lightweight standalone program containing only the inclusive
  structure-function calculation from the two codes above (no jet/analysis framework, no POWHEG
  matching). Builds two executables, `provbfh_incl` (single Higgs) and `provbfhh_incl` (di-Higgs),
  from a shared set of Fortran modules.

`releases/` holds frozen snapshots of past hepforge releases of all three programs (historical
reference — do not edit). `notes/` contains dated study logbooks (old and current) documenting
physics investigations; new studies should follow the existing pattern of a dated subdirectory with
its own `README.md`/`README.org` explaining what was done and how to reproduce it.

## Build

Each of `proVBFH/`, `proVBFHH/`, `proVBFH-inclusive/` is built independently from its own directory:

```
./configure [--with-lhapdf=/path/to/lhapdf-config] \
            [--with-hoppet=/path/to/hoppet-config] \
            [--with-fastjet=/path/to/fastjet-config]  # proVBFH/proVBFHH only
make
```

proVBFH-inclusive's `configure` differs: `--with-lhapdf=DIR` and `--with-hoppet=DIR` take the
directory containing the `*-config` binary, not the binary itself.

`configure` just writes `Makefile.inc` (external tool paths, compiler, shell); rerun it after
changing dependency locations rather than hand-editing `Makefile.inc`.

Dependencies, found via `*-config` binaries on `PATH` unless overridden:
- **hoppet** (https://github.com/hoppet-code/hoppet) — **must be v2.0.0 or newer**; the top-level
  `Makefile` in each package hard-fails the build otherwise (checked via `hoppet-config --version`).
- **LHAPDF** (http://lhapdf.hepforge.org/)
- **FastJet** (http://fastjet.fr/) — required by `proVBFH`/`proVBFHH` only, not `proVBFH-inclusive`.

`proVBFH`/`proVBFHH` build a single executable of the same name in the package root plus
`aux/combine_runs` (for merging parallel run outputs). `proVBFH-inclusive` builds both `provbfh_incl`
and `provbfhh_incl` by default (`make provbfh_incl` builds just one). `make clean` removes
object/module files; `make distclean` (proVBFH/proVBFHH) also removes the executable.
`make check` in proVBFH-inclusive builds and runs `tests/test_tensor`, a unit test of the tensor
module against its previous implementation (`tests/tensor_legacy.f90`, kept as the oracle).

Compiler defaults to gfortran (`--compiler=ifort` also supported). Object files land in
`obj-$(COMPILER)/` (proVBFH/proVBFHH) or `obj/` (proVBFH-inclusive).

## Continuous integration

`.github/workflows/ci.yml` builds all three packages on Linux and macOS against the latest releases
of hoppet, LHAPDF and FastJet, and runs short smoke tests (on push/PR, weekly, and on demand). The
logic lives in `ci/*.sh` so it can be run locally; see `ci/README.md`. The smoke tests only check
that runs complete with finite output. They are not physics validations.
`ci/tensor-checks.sh` additionally runs `make check` (proVBFH-inclusive), checks that the three
`tensor.f90` copies are identical, and compares the tensor and analytic HH matrix elements at zero
widths.

## Running

`proVBFH`/`proVBFHH` are POWHEG-style programs: run the executable from a working directory
containing a `powheg.input` and `vbfnlo.input` file (see `example/` in each package; only the inputs are
committed, and a run produces grids, `pwg-*` stage files, `xsct-*.dat` results and `.top` histogram
files). A run typically proceeds through POWHEG's staged workflow controlled by `parallelstage`/
`xgriditeration` in `powheg.input` (grid generation, then the main integration/event stage).
`aux/combine_runs` merges the outputs of multiple parallel seeds/stages.

`proVBFH-inclusive`'s executables take command-line flags directly instead of an input file, e.g.:

```
./provbfh_incl -pdf PDF4LHC15_nlo_mc -nlo -sqrts 8000 -xmur 0.5 -xmuf 0.5 -iseed 12
```

producing `xsct_(n/nn/n3)lo_seedXXXX.dat`. Valid flags are listed in `src/parameters.f90` or
`docs/provbfh-incl-doc.pdf`.

Documentation for each package is written in LaTeX under `docs/*.tex`, with the built PDF checked in
alongside it.

## Code architecture (proVBFH / proVBFHH)

Source is split by role, reflected both in directory layout and in the Makefile's object-file
groupings (`INCL`, `EXCL`, `PWHG`, `HJJJ`, `VBFNLO`):

- **`src/inclusive/`** — the structure-function matrix elements and phase space used for the
  inclusive/NNLO+ part of the calculation (`incl_vbfh.f90`, `matrix_element.f90`,
  `incl_parameters.f90`, `phase_space.f`, `tensor.f90`, `nonfact.f90` for non-factorisable
  corrections). This is the same physics content that `proVBFH-inclusive` reuses as a standalone
  program.
- **`src/exclusive/`** — glues the inclusive structure functions into the POWHEG BOX (`excl_vbfh.f90`,
  `Born.f`, `real_vbfnlo.f`, `sigreal.f`, `pwhg_init.f`, ...), including a bundled/patched copy of
  VBFNLO's Hjjj matrix elements under `exclusive/vbfnlo-files`.
- **`src/powheg-files/`** — the generic POWHEG BOX framework files (phase space generation, grid
  handling, LHE output, `powheginput.f`, etc.), largely shared boilerplate across POWHEG processes;
  `powheg-files/hjjj-files` holds the Hjjj-specific virtual/real matrix element pieces.
- **`analysis/`** — user analysis routines (histogramming, cuts) selected via the `ANALYSIS_FILE`
  variable in the Makefile (only one analysis is linked in at a time — swap by editing that variable).

`src/inclusive/tensor.f90` (proVBFH, proVBFHH) and `src/tensor.f90` (proVBFH-inclusive) are the same
file; keep the three copies identical. Four-vectors go into tensors with `InitFourVector(t, p, up)`
(contravariant `p`, lowered with the metric if `up` is false); storing `p` directly in a tensor with
lower indices describes the parity-flipped vector, which broke the F3 terms of the HH tensor matrix
element before 03adc62.

Module/compile-order dependencies are declared explicitly at the bottom of the Makefile (e.g.
`incl_vbfh.o: matrix_element.o incl_parameters.o phase_space.o ...`); when adding new `.f90` modules,
add the corresponding dependency line or the build can fail/use stale `.mod` files.

`Makefile.inc` (generated by `configure`) is included by `Makefile` and is the only place external
tool paths/compiler choice should be changed.

## Code architecture (proVBFH-inclusive)

Much flatter: all sources live directly in `src/`. `provbfh_incl.f90` and `provbfhh_incl.f90` are the
two main programs, sharing modules `parameters.f90`, `integration.f`, `io_utils.f90`, `lcl_dec.f90`,
`tensor.f90`; single-Higgs adds `matrix_element.f90`/`phase_space.f`, di-Higgs adds
`matrix_element_dihiggs.f90`/`phase_space_dihiggs.f`/`ME_expressions.f`. As in the other packages,
module dependencies are listed explicitly in the Makefile and must be kept in sync when the module
graph changes.

Several `*_new.f90`/`*_new.f` and `*_save.f90` files alongside their originals in this directory are
work-in-progress/backup variants (e.g. `parameters_new.f90` vs `parameters.f90`) rather than files
currently wired into the Makefile — check `Makefile`'s `MODULES`/`MAIN` lists before assuming a `.f90`
file is part of the active build.
