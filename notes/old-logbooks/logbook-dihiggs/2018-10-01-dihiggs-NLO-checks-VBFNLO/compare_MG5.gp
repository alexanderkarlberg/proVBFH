set term pdfcairo enhanced color lw 2 size 4.0in,3.0in
set macros
set datafile fortran
set grid

fix=' | sed -e "s/D/E/g" | tail -n +7'
pbtofb=1

fnLO_pr='proVBFH/LO-total.dat '
fnLO_mg='MG5/LO-hists.top '
fnNLO_pr='proVBFH/NLO-total.dat '
fnNLO_mg='MG5/NLO-hists.top '

ptja_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptj1)  <(mergeidx.pl -td -f '.fnLO_mg.'ptj1.vbf '.fix.')"'
ptjb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptj2)  <(mergeidx.pl -td -f '.fnLO_mg.'ptj2.vbf '.fix.')"'
ptHd_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptH_h) <(mergeidx.pl -td -f '.fnLO_mg.'ptH_h.vbf'.fix.')"'
ptHs_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptH_s) <(mergeidx.pl -td -f '.fnLO_mg.'ptH_s.vbf'.fix.')"'
ptHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'ptHH)  <(mergeidx.pl -td -f '.fnLO_mg.'ptHH.vbf '.fix.')"'

ptja_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptj1)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptj1.vbf '.fix.')"'
ptjb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptj2)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptj2.vbf '.fix.')"'
ptHd_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptH_h) <(mergeidx.pl -td -f '.fnNLO_mg.'ptH_h.vbf'.fix.')"'
ptHs_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptH_s) <(mergeidx.pl -td -f '.fnNLO_mg.'ptH_s.vbf'.fix.')"'
ptHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'ptHH)  <(mergeidx.pl -td -f '.fnNLO_mg.'ptHH.vbf '.fix.')"'

yja_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yj1)  <(mergeidx.pl -td -f '.fnLO_mg.'yj1.vbf '.fix.')"'
yjb_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yj2)  <(mergeidx.pl -td -f '.fnLO_mg.'yj2.vbf '.fix.')"'
yHd_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yH_h) <(mergeidx.pl -td -f '.fnLO_mg.'yH_h.vbf'.fix.')"'
yHs_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yH_s) <(mergeidx.pl -td -f '.fnLO_mg.'yH_s.vbf'.fix.')"'
yHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'yHH)  <(mergeidx.pl -td -f '.fnLO_mg.'yHH.vbf '.fix.')"'

yja_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yj1)  <(mergeidx.pl -td -f '.fnNLO_mg.'yj1.vbf '.fix.')"'
yjb_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yj2)  <(mergeidx.pl -td -f '.fnNLO_mg.'yj2.vbf '.fix.')"'
yHd_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yH_h) <(mergeidx.pl -td -f '.fnNLO_mg.'yH_h.vbf'.fix.')"'
yHs_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yH_s) <(mergeidx.pl -td -f '.fnNLO_mg.'yH_s.vbf'.fix.')"'
yHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'yHH)  <(mergeidx.pl -td -f '.fnNLO_mg.'yHH.vbf '.fix.')"'

mHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'mHH) <(mergeidx.pl -td -f '.fnLO_mg.'mHH.vbf'.fix.')"'
phiHH_LO = '< exec bash -c "paste <(mergeidx.pl -f '.fnLO_pr.'phi.H1.H2) <(mergeidx.pl -td -f '.fnLO_mg.'DphiH1H2.vbf'.fix.')"'

mHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'mHH) <(mergeidx.pl -td -f '.fnNLO_mg.'MHH.vbf'.fix.')"'
phiHH_NLO = '< exec bash -c "paste <(mergeidx.pl -f '.fnNLO_pr.'phi.H1.H2) <(mergeidx.pl -td -f '.fnNLO_mg.'DphiH1H2.vbf'.fix.')"'


print ptHH_LO

bpta=20
bptb=10
bpt=30
byH=0.1
by=0.5
bphi=0.31415

set output 'comparison-MG5-LO.pdf'
set xlabel 'p_{t,HH}'
plot ptHH_LO u 5:($6/bpt) w l t 'MadGraph',\
     ptHH_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

print yHH_LO
set xlabel 'y_{HH}'
plot yHH_LO u 5:($6/byH) w l t 'MadGraph',\
     yHH_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH',\

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_LO u 5:($6/bphi) w l t 'MadGraph',\
     phiHH_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,j_1}'
plot ptja_LO u 5:($6/bpta) w l t 'MadGraph',\
     ptja_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_1}'
plot yja_LO u 5:($6/by) w l t 'MadGraph',\
     yja_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,j_2}'
plot ptjb_LO u 5:($6/bptb) w l t 'MadGraph',\
     ptjb_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{j_2}'
plot yjb_LO u 5:($6/by) w l t 'MadGraph',\
     yjb_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_LO u 5:($6/bpt) w l t 'MadGraph',\
     ptHd_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{hard}}'
plot yHd_LO u 5:($6/byH) w l t 'MadGraph',\
     yHd_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_LO u 5:($6/bpt) w l t 'MadGraph',\
     ptHs_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set xlabel 'y_{H_{soft}}'
plot yHs_LO u 5:($6/byH) w l t 'MadGraph',\
     yHs_LO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

set output 'comparison-MG5-ratio-LO.pdf'
set yrange [0.97:1.03]
set xlabel 'p_{t,HH}'
plot ptHH_LO u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{HH}'
plot yHH_LO u (0.5*($1+$2)):(pbtofb*$3*byH/$6) w l t 'proVBFH/MadGraph'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_LO u (0.5*($1+$2)):(pbtofb*$3*bphi/$6) w l t 'proVBFH/MadGraph'

set xlabel 'p_{t,j_1}'
plot ptja_LO u (0.5*($1+$2)):(pbtofb*$3*bpta/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{j_1}'
plot yja_LO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel 'p_{t,j_2}'
plot ptjb_LO u (0.5*($1+$2)):(pbtofb*$3*bptb/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{j_2}'
plot yjb_LO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_LO u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{H_{hard}}'
plot yHd_LO u (0.5*($1+$2)):(pbtofb*$3*byH/$6)   w l t 'proVBFH/MadGraph'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_LO u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MadGraph'

set xlabel 'y_{H_{soft}}'
plot yHs_LO u (0.5*($1+$2)):(pbtofb*$3*byH/$6)   w l t 'proVBFH/MadGraph'


# set autoscale y
# set output 'comparison-MG5-NLO.pdf'
# set xlabel 'p_{t,HH}'
# plot ptHH_NLO u 5:($6/bpt) w l t 'MadGraph',\
#      ptHH_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# print yHH_NLO
# set xlabel 'y_{HH}'
# plot yHH_NLO u 5:($6/byH) w l t 'MadGraph',\
#      yHH_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH',\

# set xlabel '{/Symbol D}{/Symbol f}_{HH}'
# plot phiHH_NLO u 5:($6/bphi) w l t 'MadGraph',\
#      phiHH_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'p_{t,j_1}'
# plot ptja_NLO u 5:($6/bpta) w l t 'MadGraph',\
#      ptja_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'y_{j_1}'
# plot yja_NLO u 5:($6/by) w l t 'MadGraph',\
#      yja_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'p_{t,j_2}'
# plot ptjb_NLO u 5:($6/bptb) w l t 'MadGraph',\
#      ptjb_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'y_{j_2}'
# plot yjb_NLO u 5:($6/by) w l t 'MadGraph',\
#      yjb_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'p_{t,H_{hard}}'
# plot ptHd_NLO u 5:($6/bpt) w l t 'MadGraph',\
#      ptHd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'y_{H_{hard}}'
# plot yHd_NLO u 5:($6/byH) w l t 'MadGraph',\
#      yHd_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'p_{t,H_{soft}}'
# plot ptHs_NLO u 5:($6/bpt) w l t 'MadGraph',\
#      ptHs_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set xlabel 'y_{H_{soft}}'
# plot yHs_NLO u 5:($6/byH) w l t 'MadGraph',\
#      yHs_NLO u (0.5*($1+$2)):(pbtofb*$3) w l dt 2 t 'proVBFH'

# set output 'comparison-MG5-ratio-NLO.pdf'
# set yrange [0.92:1.08]
# set xlabel 'p_{t,HH}'
# plot ptHH_NLO u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MadGraph'

# set xlabel 'y_{HH}'
# plot yHH_NLO u (0.5*($1+$2)):(pbtofb*$3*byHH/$6) w l t 'proVBFH/MadGraph'

# set xlabel '{/Symbol D}{/Symbol f}_{HH}'
# plot phiHH_NLO u (0.5*($1+$2)):(pbtofb*$3*bphi/$6) w l t 'proVBFH/MadGraph'

# set xlabel 'p_{t,j_1}'
# plot ptja_NLO u (0.5*($1+$2)):(pbtofb*$3*bpta/$6) w l t 'proVBFH/MadGraph'

# set xlabel 'y_{j_1}'
# plot yja_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

# set xlabel 'p_{t,j_2}'
# plot ptjb_NLO u (0.5*($1+$2)):(pbtofb*$3*bptb/$6) w l t 'proVBFH/MadGraph'

# set xlabel 'y_{j_2}'
# plot yjb_NLO u (0.5*($1+$2)):(pbtofb*$3*by/$6)   w l t 'proVBFH/MadGraph'

# set xlabel 'p_{t,H_{hard}}'
# plot ptHd_NLO u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MadGraph'

# set xlabel 'y_{H_{hard}}'
# plot yHd_NLO u (0.5*($1+$2)):(pbtofb*$3*byH/$6)   w l t 'proVBFH/MadGraph'

# set xlabel 'p_{t,H_{soft}}'
# plot ptHs_NLO u (0.5*($1+$2)):(pbtofb*$3*bpt/$6) w l t 'proVBFH/MadGraph'

# set xlabel 'y_{H_{soft}}'
# plot yHs_NLO u (0.5*($1+$2)):(pbtofb*$3*byH/$6)   w l t 'proVBFH/MadGraph'

