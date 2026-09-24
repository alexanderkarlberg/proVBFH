# F3 check: index-correct calculation vs the two HH matrix elements

Evaluates the squared VBF HH matrix element contracted with the F3
parts of the hadronic tensors at one physical phase-space point
(momentum conservation, on-shell Higgs bosons), with explicit
upper/lower index handling (`lowerit=True`), and with the convention
of the tensor code before 03adc62 (four-vector components stored
unlowered in tensors with lower indices, `lowerit=False`).

    python3 check.py            # correct vs old-tensor convention
    python3 coeffs.py 1         # F1F1 and F3F3 coefficients AA..CC, writes point.txt
    gfortran -ffixed-line-length-132 -c ../../../../proVBFH-inclusive/src/ME_expressions.f -o ME.o
    gfortran drv.f90 ME.o -o drv && ./drv   # analytic F1F1, F3F3 at point.txt
    python3 width.py            # F1F3, F2F3 terms with and without widths

With seed 1 the analytic F3F3 coefficients agree with the correct
calculation to all printed digits; the old-tensor convention agrees
for AA only.
