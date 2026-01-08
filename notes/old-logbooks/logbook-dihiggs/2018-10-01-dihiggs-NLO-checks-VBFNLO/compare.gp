set term pdfcairo enhanced color lw 2 size 4.0in,3.0in
set macros
set datafile fortran
set grid

pbtofb=1000

fnLO_pr='proVBFH/LO-total.dat '
fnLO_vb(n)='VBFNLO/LO/hist.'.n.'.dat'
fnNLO_pr='proVBFH/NLO-total.dat '
fnNLO_vb(n)='VBFNLO/NLO/hist.'.n.'.dat'

ptHa_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptH1)  '.fnLO_vb('16').'"'
ptHb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptH2)  '.fnLO_vb('17').'"'
ptHd_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptH_h) '.fnLO_vb('18').'"'
ptHs_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptH_s) '.fnLO_vb('19').'"'
ptHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptHH)  '.fnLO_vb('15').'"'

ptHa_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptH1)  '.fnNLO_vb('16').'"'
ptHb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptH2)  '.fnNLO_vb('17').'"'
ptHd_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptH_h) '.fnNLO_vb('18').'"'
ptHs_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptH_s) '.fnNLO_vb('19').'"'
ptHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptHH)  '.fnNLO_vb('15').'"'

yHa_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yH1)  '.fnLO_vb('21').'"'
yHb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yH2)  '.fnLO_vb('22').'"'
yHd_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yH_h) '.fnLO_vb('23').'"'
yHs_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yH_s) '.fnLO_vb('24').'"'
yHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yHH)  '.fnLO_vb('20').'"'

yHa_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yH1)  '.fnNLO_vb('21').'"'
yHb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yH2)  '.fnNLO_vb('22').'"'
yHd_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yH_h) '.fnNLO_vb('23').'"'
yHs_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yH_s) '.fnNLO_vb('24').'"'
yHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yHH)  '.fnNLO_vb('20').'"'

mHa_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'mH1) '.fnLO_vb('13').'"'
mHb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'mH2) '.fnLO_vb('14').'"'
mHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'mHH) '.fnLO_vb('12').'"'
phiHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'phi.H1.H2) '.fnLO_vb('25').'"'

mHa_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'mH1) '.fnNLO_vb('13').'"'
mHb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'mH2) '.fnNLO_vb('14').'"'
mHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'mHH) '.fnNLO_vb('12').'"'
phiHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'phi.H1.H2) '.fnNLO_vb('25').'"'


print ptHH_LO
set output 'comparison-vbfnlo-LO.pdf'
set xlabel 'p_{t,HH}'
plot ptHH_LO u 5:6 w l t 'VBFNLO',\
     ptHH_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

print yHH_LO
set xlabel 'y_{HH}'
plot yHH_LO u 5:6 w l t 'VBFNLO',\
     yHH_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH',\

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_LO u 5:6 w l t 'VBFNLO',\
     phiHH_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_1}'
plot ptHa_LO u 5:6 w l t 'VBFNLO',\
     ptHa_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_1}'
plot yHa_LO u 5:6 w l t 'VBFNLO',\
     yHa_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_2}'
plot ptHb_LO u 5:6 w l t 'VBFNLO',\
     ptHb_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_2}'
plot yHb_LO u 5:6 w l t 'VBFNLO',\
     yHb_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_LO u 5:6 w l t 'VBFNLO',\
     ptHd_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{hard}}'
plot yHd_LO u 5:6 w l t 'VBFNLO',\
     yHd_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_LO u 5:6 w l t 'VBFNLO',\
     ptHs_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{soft}}'
plot yHs_LO u 5:6 w l t 'VBFNLO',\
     yHs_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set output 'comparison-vbfnlo-ratio-LO.pdf'
set yrange [0.97:1.03]
set xlabel 'p_{t,HH}'
plot ptHH_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{HH}'
plot yHH_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_1}'
plot ptHa_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_1}'
plot yHa_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_2}'
plot ptHb_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_2}'
plot yHb_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_{hard}}'
plot yHd_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_{soft}}'
plot yHs_LO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'


set autoscale y
set output 'comparison-vbfnlo-NLO.pdf'
set xlabel 'p_{t,HH}'
plot ptHH_NLO u 5:6 w l t 'VBFNLO',\
     ptHH_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

print yHH_NLO
set xlabel 'y_{HH}'
plot yHH_NLO u 5:6 w l t 'VBFNLO',\
     yHH_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH',\

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_NLO u 5:6 w l t 'VBFNLO',\
     phiHH_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_1}'
plot ptHa_NLO u 5:6 w l t 'VBFNLO',\
     ptHa_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_1}'
plot yHa_NLO u 5:6 w l t 'VBFNLO',\
     yHa_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_2}'
plot ptHb_NLO u 5:6 w l t 'VBFNLO',\
     ptHb_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_2}'
plot yHb_NLO u 5:6 w l t 'VBFNLO',\
     yHb_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_NLO u 5:6 w l t 'VBFNLO',\
     ptHd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{hard}}'
plot yHd_NLO u 5:6 w l t 'VBFNLO',\
     yHd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_NLO u 5:6 w l t 'VBFNLO',\
     ptHs_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{soft}}'
plot yHs_NLO u 5:6 w l t 'VBFNLO',\
     yHs_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set output 'comparison-vbfnlo-ratio-NLO.pdf'
set yrange [0.92:1.08]
set xlabel 'p_{t,HH}'
plot ptHH_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{HH}'
plot yHH_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_1}'
plot ptHa_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_1}'
plot yHa_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_2}'
plot ptHb_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_2}'
plot yHb_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_{hard}}'
plot yHd_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

set xlabel 'y_{H_{soft}}'
plot yHs_NLO u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'proVBFH/VBFNLO'

