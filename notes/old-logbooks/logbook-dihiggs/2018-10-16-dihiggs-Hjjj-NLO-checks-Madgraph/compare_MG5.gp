set term pdfcairo enhanced color lw 2 size 4.0in,3.0in
set macros
set datafile fortran
set grid

fix=' | sed -e "s/D/E/g" | tail -n +7'
pbtofb=1

fnLO_pr='proVBFH/Hjjj-LO.top '
fnLO_mg='MG5/Hjjj-LO.top '
fnNLO_pr='proVBFH/Hjjj-NLO.top '
fnNLO_mg='MG5/Hjjj-NLO.top '

ptja_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptj1)  <(mergeidx.pl -td -f '.fnLO_mg.'ptj1.vbf '.fix.')"'
ptjb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptj2)  <(mergeidx.pl -td -f '.fnLO_mg.'ptj2.vbf '.fix.')"'
ptjc_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptj3)  <(mergeidx.pl -td -f '.fnLO_mg.'ptj3.vbf '.fix.')"'

ptja_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptj1)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptj1.vbf '.fix.')"'
ptjb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptj2)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptj2.vbf '.fix.')"'
ptjc_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptj3)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptj3.vbf '.fix.')"'
ptjd_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptj4)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptj4.vbf '.fix.')"'

yja_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yj1)  <(mergeidx.pl -td -f '.fnLO_mg.'yj1.vbf '.fix.')"'
yjb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yj2)  <(mergeidx.pl -td -f '.fnLO_mg.'yj2.vbf '.fix.')"'
yjc_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yj3)  <(mergeidx.pl -td -f '.fnLO_mg.'yj3.vbf '.fix.')"'
yjcs_LO= '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'y.j3)  <(mergeidx.pl -td -f '.fnLO_mg.'y.j3.vbf '.fix.')"'
yjac_LO= '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'rapj1j3)  <(mergeidx.pl -td -f '.fnLO_mg.'rapj1j3.vbf '.fix.')"'
yjbc_LO= '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'rapj2j3)  <(mergeidx.pl -td -f '.fnLO_mg.'rapj2j3.vbf '.fix.')"'

yja_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yj1)  <(mergeidx.pl -td -f '.fnNLO_mg.'yj1.vbf '.fix.')"'
yjb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yj2)  <(mergeidx.pl -td -f '.fnNLO_mg.'yj2.vbf '.fix.')"'
yjc_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yj3)  <(mergeidx.pl -td -f '.fnNLO_mg.'yj3.vbf '.fix.')"'
yjd_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yj4)  <(mergeidx.pl -td -f '.fnNLO_mg.'yj4.vbf '.fix.')"'
yjcs_NLO= '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'y.j3)  <(mergeidx.pl -td -f '.fnNLO_mg.'y.j3.vbf '.fix.')"'
yjac_NLO= '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'rapj1j3)  <(mergeidx.pl -td -f '.fnNLO_mg.'rapj1j3.vbf '.fix.')"'
yjbc_NLO= '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'rapj2j3)  <(mergeidx.pl -td -f '.fnNLO_mg.'rapj2j3.vbf '.fix.')"'
yjad_NLO= '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'rapj1j4)  <(mergeidx.pl -td -f '.fnNLO_mg.'rapj1j4.vbf '.fix.')"'
yjbd_NLO= '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'rapj2j4)  <(mergeidx.pl -td -f '.fnNLO_mg.'rapj2j4.vbf '.fix.')"'
yjcd_NLO= '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'rapj3j4)  <(mergeidx.pl -td -f '.fnNLO_mg.'rapj3j4.vbf '.fix.')"'

print ptja_LO

bpta=20
bptb=10
bptc=5
bpt=30
byH=0.1
by=0.5
bphi=0.31415

set output 'comparison-MG5-LO.pdf'

set xlabel 'p_{t,j_3}'
plot ptjc_LO u 5:($6/bptc) w l t 'MadGraph',\
     ptjc_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_3}'
plot yjc_LO u 5:($6/by) w l t 'MadGraph',\
     yjc_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y^*_{j_3}'
plot yjcs_LO u 5:($6/by) w l t 'MadGraph',\
     yjcs_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_1} - y_{j_3}|'
plot yjac_LO u 5:($6/by) w l t 'MadGraph',\
     yjac_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_2} - y_{j_3}|'
plot yjbc_LO u 5:($6/by) w l t 'MadGraph',\
     yjbc_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set output 'comparison-MG5-ratio-LO.pdf'
set yrange [0.92:1.08]

set xlabel 'p_{t,j_3}'
plot ptjc_LO u (0.5*($1+$2)):(pbtofb*$3*bptc/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{j_3}'
plot yjc_LO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel 'y^*_{j_3}'
plot yjcs_LO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_1} - y_{j_3}|'
plot yjac_LO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_2} - y_{j_3}|'
plot yjbc_LO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'


reset


set output 'comparison-MG5-NLO.pdf'

set xlabel 'p_{t,j_3}'
plot ptjc_NLO u 5:($6/bptc) w l t 'MadGraph',\
     ptjc_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_3}'
plot yjc_NLO u 5:($6/by) w l t 'MadGraph',\
     yjc_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y^*_{j_3}'
plot yjcs_NLO u 5:($6/by) w l t 'MadGraph',\
     yjcs_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_1} - y_{j_3}|'
plot yjac_NLO u 5:($6/by) w l t 'MadGraph',\
     yjac_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_2} - y_{j_3}|'
plot yjbc_NLO u 5:($6/by) w l t 'MadGraph',\
     yjbc_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,j_4}'
plot ptjd_NLO u 5:($6/bptc) w l t 'MadGraph',\
     ptjd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_4}'
plot yjd_NLO u 5:($6/by) w l t 'MadGraph',\
     yjd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_1} - y_{j_4}|'
plot yjad_NLO u 5:($6/by) w l t 'MadGraph',\
     yjad_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_2} - y_{j_4}|'
plot yjbd_NLO u 5:($6/by) w l t 'MadGraph',\
     yjbd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel '|y_{j_3} - y_{j_4}|'
plot yjcd_NLO u 5:($6/by) w l t 'MadGraph',\
     yjcd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set output 'comparison-MG5-ratio-NLO.pdf'
set yrange [0.7:1.3]

set xlabel 'p_{t,j_3}'
plot ptjc_NLO u (0.5*($1+$2)):(pbtofb*$3*bptc/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{j_3}'
plot yjc_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel 'y^*_{j_3}'
plot yjcs_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_1} - y_{j_3}|'
plot yjac_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_2} - y_{j_3}|'
plot yjbc_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel 'p_{t,j_4}'
plot ptjd_NLO u (0.5*($1+$2)):(pbtofb*$3*bptc/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{j_4}'
plot yjd_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_1} - y_{j_4}|'
plot yjad_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_2} - y_{j_4}|'
plot yjbd_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel '|y_{j_3} - y_{j_4}|'
plot yjcd_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

