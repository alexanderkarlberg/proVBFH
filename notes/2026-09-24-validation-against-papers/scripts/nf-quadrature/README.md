# Quadrature study for the NF azimuthal integrals

`args.dat.gz`: 9000 argument sets (tag, MV, MVH2, p1x, p2x, p2y, p3x,
p3y) of `box_1loop_new` (B1), `box_2loop` (B2) and `tri_2loop` (T2),
dumped from the first iteration of HH and H NF runs (2005.11334
setup). The programs evaluate the integrands of
`proVBFHH/src/inclusive/nonfact_expressions.f` at these points:

    gunzip -k args.dat.gz
    E=../../../../proVBFHH/src/inclusive/nonfact_expressions.f
    git show f37f0da~1:proVBFHH/src/inclusive/nonfact_expressions.f > ne_sp.f   # single precision
    sed 's/module nonfact_expressions/module nonfact_expressions_dp/; s/complex\*8/complex*16/' ne_sp.f > ne_dp.f
    gfortran -O2 -ffixed-line-length-132 -c ne_sp.f ne_dp.f && gfortran -O2 -c integ.f90
    gfortran -O2 ne_sp.o ne_dp.o integ.o quad.f90 -o quad && ./quad    # RK4 vs trapezoid vs Gauss-Legendre
    gfortran -O2 ne_sp.o ne_dp.o integ.o gk.f90 -o gk && ./gk          # adaptive Gauss-Kronrod
    gfortran -O2 ne_sp.o ne_dp.o integ.o worst.f90 -o worst && ./worst # worst t022 point, error quantiles
    gfortran -O2 ne_sp.o ne_dp.o integ.o prec.f90 -o prec && ./prec    # single vs double precision
    gfortran -O2 -ffixed-line-length-132 $E pair.f90 -o pair && ./pair # box_integrands vs b01, b022
