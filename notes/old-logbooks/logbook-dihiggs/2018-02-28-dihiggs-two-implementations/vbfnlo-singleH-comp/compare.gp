set term pdfcairo enhanced color lw 2 size 4.0in,3.0in
set macros
set datafile fortran
set grid

pbtofb=1000

ptH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top   ptH) histograms.dir/LO/hist.15.dat"'
yH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top   yH) histograms.dir/LO/hist.20.dat"'
mH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top   mH) histograms.dir/LO/hist.12.dat"'

set output 'comparison-vbfnlo.pdf'
set xlabel 'p_{t,H}'
plot ptH u 5:6 w l t 'VBFNLO',\
     ptH u (0.5*($1+$2)):(pbtofb*$3) w l t 'proVBFH'

set xlabel 'y_{H}'
plot yH u 5:6 w l t 'VBFNLO',\
     yH u (0.5*($1+$2)):(pbtofb*$3) w l t 'proVBFH'

set output 'comparison-vbfnlo-ratio.pdf'
set yrange [0.92:1.08]
set xlabel 'p_{t,H}'
plot ptH u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H}'
plot yH u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

