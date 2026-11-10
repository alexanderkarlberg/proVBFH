# Exact-check references

Frozen outputs of short, fixed proVBFH-cs runs, compared byte for byte
by `ci/physics-checks.py exact` (see its docstring). Regenerate with
`ci/physics-checks.py generate exact` after a deliberate change of the
code or of the dependencies or compiler.

- Generated: 2026-11-10
- Code: commit e1c3dac
- gfortran: GNU Fortran (Ubuntu 15.2.0-16ubuntu1) 15.2.0
- hoppet: 2.3.1
- LHAPDF: 6.5.6
- PDF: the hoppet benchmark (ci/benchmark-pdf.sh)
- Runs: qcd_order 2 and 3, iseed 1, ncall1 = ncall2 = 2000, itmx1 = itmx2 = 1,
  13 TeV, scale of 1506.02660 (cs_scales 3), default card otherwise
- Compared: xsct-*.dat and pwg-*.top of each run (not the grids)
