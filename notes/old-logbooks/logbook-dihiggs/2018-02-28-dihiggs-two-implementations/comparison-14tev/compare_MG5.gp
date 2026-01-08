set term pdfcairo enhanced color lw 2 size 4.0in,3.0in
set macros
set datafile fortran
set grid

pbtofb=1
fix=' | sed -e "s/D/E/g" | tail -n +7'

ptja = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top ptj1)  <(mergeidx.pl -td -f mg5-LO.top ptj1'.fix.')"'
ptjb = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top ptj2)  <(mergeidx.pl -td -f mg5-LO.top ptj2'.fix.')"'
ptHd = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top ptH_h) <(mergeidx.pl -td -f mg5-LO.top ptH_h'.fix.')"'
ptHs = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top ptH_s) <(mergeidx.pl -td -f mg5-LO.top ptH_s'.fix.')"'
ptHH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top ptHH)  <(mergeidx.pl -td -f mg5-LO.top ptHH'.fix.')"'

yja = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top yj1)    <(mergeidx.pl -td -f mg5-LO.top yj1 '.fix.')"'
yjb = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top yj2)    <(mergeidx.pl -td -f mg5-LO.top yj2 '.fix.')"'
yHd = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top yH_h)   <(mergeidx.pl -td -f mg5-LO.top yH_h'.fix.')"'
yHs = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top yH_s)   <(mergeidx.pl -td -f mg5-LO.top yH_s'.fix.')"'
yHH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top yHH)    <(mergeidx.pl -td -f mg5-LO.top yHH '.fix.')"'

mHa = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top mH1)    <(mergeidx.pl -td -f mg5-LO.top mH1 '.fix.')"'
mHb = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top mH2)    <(mergeidx.pl -td -f mg5-LO.top mH2 '.fix.')"'
mHH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top mHH)    <(mergeidx.pl -td -f mg5-LO.top mHH '.fix.')"'

phiHH = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO.top phi.H1.H2) <(mergeidx.pl -td -f mg5-LO.top DphiH1H2'.fix.')"'

bpta=20
bptb=10
bpt=30
byH=0.1
by=0.5
bphi=0.31415
set output 'comparison-madgraph.pdf'
set xlabel 'p_{t,HH}'
plot ptHH u 5:($6/bpt) w l t 'MADGRAPH',\
     ptHH u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{HH}'
plot yHH u 5:($6/byH) w l t 'MADGRAPH',\
     yHH u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH u 5:($6/bphi) w l t 'MADGRAPH',\
     phiHH u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,j_1}'
plot ptja u 5:($6/bpta) w l t 'MADGRAPH',\
     ptja u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_1}'
plot yja u 5:($6/by) w l t 'MADGRAPH',\
     yja u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,j_2}'
plot ptjb u 5:($6/bptb) w l t 'MADGRAPH',\
     ptjb u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_2}'
plot yjb u 5:($6/by) w l t 'MADGRAPH',\
     yjb u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{hard}}'
plot ptHd u 5:($6/bpt) w l t 'MADGRAPH',\
     ptHd u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{hard}}'
plot yHd u 5:($6/byH) w l t 'MADGRAPH',\
     yHd u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{soft}}'
plot ptHs u 5:($6/bpt) w l t 'MADGRAPH',\
     ptHs u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{soft}}'
plot yHs u 5:($6/byH) w l t 'MADGRAPH',\
     yHs u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set output 'comparison-madgraph-ratio.pdf'
set yrange [0.92:1.08]
set xlabel 'p_{t,HH}'
plot ptHH u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MADGRAPH'

set xlabel 'y_{HH}'
plot yHH u (0.5*($1+$2)):(pbtofb*$3*byH/$6) w l t 'proVBFH/MADGRAPH'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH u (0.5*($1+$2)):(pbtofb*$3*bphi/$6) w l t 'proVBFH/MADGRAPH'

set xlabel 'p_{t,j_1}'
plot ptja u (0.5*($1+$2)):(pbtofb*$3*bpta/$6) w l t 'proVBFH/MADGRAPH'

set xlabel 'y_{j_1}'
plot yja u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MADGRAPH'

set xlabel 'p_{t,j_2}'
plot ptjb u (0.5*($1+$2)):(pbtofb*$3*bptb/$6) w l t 'proVBFH/MADGRAPH'

set xlabel 'y_{j_2}'
plot yjb u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MADGRAPH'

set xlabel 'p_{t,H_{hard}}'
plot ptHd u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MADGRAPH'

set xlabel 'y_{H_{hard}}'
plot yHd u (0.5*($1+$2)):(pbtofb*$3*byH/$6)   w l t 'proVBFH/MADGRAPH'

set xlabel 'p_{t,H_{soft}}'
plot ptHs u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MADGRAPH'

set xlabel 'y_{H_{soft}}'
plot yHs u (0.5*($1+$2)):(pbtofb*$3*byH/$6)   w l t 'proVBFH/MADGRAPH'

