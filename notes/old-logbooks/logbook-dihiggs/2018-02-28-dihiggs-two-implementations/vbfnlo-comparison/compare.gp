set term pdfcairo enhanced color lw 2 size 4.0in,3.0in
set macros
set datafile fortran
set grid

pbtofb=1000

ptHa_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   ptH1) histograms.dir/LO/hist.16.dat"'
ptHa_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top ptH1) histograms.dir/LO/hist.16.dat"'
ptHb_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   ptH2) histograms.dir/LO/hist.17.dat"'
ptHb_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top ptH2) histograms.dir/LO/hist.17.dat"'
ptHd_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   ptH_h) histograms.dir/LO/hist.18.dat"'
ptHd_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top ptH_h) histograms.dir/LO/hist.18.dat"'
ptHs_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   ptH_s) histograms.dir/LO/hist.19.dat"'
ptHs_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top ptH_s) histograms.dir/LO/hist.19.dat"'
ptHH_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   ptHH) histograms.dir/LO/hist.15.dat"'
ptHH_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top ptHH) histograms.dir/LO/hist.15.dat"'

yHa_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   yH1) histograms.dir/LO/hist.21.dat"'
yHa_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top yH1) histograms.dir/LO/hist.21.dat"'
yHb_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   yH2) histograms.dir/LO/hist.22.dat"'
yHb_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top yH2) histograms.dir/LO/hist.22.dat"'
yHd_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   yH_h) histograms.dir/LO/hist.23.dat"'
yHd_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top yH_h) histograms.dir/LO/hist.23.dat"'
yHs_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   yH_s) histograms.dir/LO/hist.24.dat"'
yHs_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top yH_s) histograms.dir/LO/hist.24.dat"'
yHH_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   yHH) histograms.dir/LO/hist.20.dat"'
yHH_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top yHH) histograms.dir/LO/hist.20.dat"'

mHa_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   mH1) histograms.dir/LO/hist.13.dat"'
mHa_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top mH1) histograms.dir/LO/hist.13.dat"'
mHb_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   mH2) histograms.dir/LO/hist.14.dat"'
mHb_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top mH2) histograms.dir/LO/hist.14.dat"'
mHH_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   mHH) histograms.dir/LO/hist.12.dat"'
mHH_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top mHH) histograms.dir/LO/hist.12.dat"'

phiHH_ten = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-tensor.top   phi.H1.H2) histograms.dir/LO/hist.25.dat"'
phiHH_ana = '< exec bash -c "paste <(mergeidx.pl -f pwg-LO-analytic.top phi.H1.H2) histograms.dir/LO/hist.25.dat"'

set output 'comparison-vbfnlo.pdf'
set xlabel 'p_{t,HH}'
plot ptHH_ana u 5:6 w l t 'VBFNLO',\
     ptHH_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     ptHH_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'y_{HH}'
plot yHH_ana u 5:6 w l t 'VBFNLO',\
     yHH_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     yHH_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'm_{HH}'
plot mHH_ana u 5:6 w l t 'VBFNLO',\
     mHH_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     mHH_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_ana u 5:6 w l t 'VBFNLO',\
     phiHH_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     phiHH_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'p_{t,H_1}'
plot ptHa_ana u 5:6 w l t 'VBFNLO',\
     ptHa_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     ptHa_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'y_{H_1}'
plot yHa_ana u 5:6 w l t 'VBFNLO',\
     yHa_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     yHa_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'm_{H_1}'
plot mHa_ana u 5:6 w l t 'VBFNLO',\
     mHa_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     mHa_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'p_{t,H_2}'
plot ptHb_ana u 5:6 w l t 'VBFNLO',\
     ptHb_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     ptHb_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'y_{H_2}'
plot yHb_ana u 5:6 w l t 'VBFNLO',\
     yHb_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     yHb_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'm_{H_2}'
plot mHb_ana u 5:6 w l t 'VBFNLO',\
     mHb_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     mHb_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_ana u 5:6 w l t 'VBFNLO',\
     ptHd_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     ptHd_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'y_{H_{hard}}'
plot yHd_ana u 5:6 w l t 'VBFNLO',\
     yHd_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     yHd_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_ana u 5:6 w l t 'VBFNLO',\
     ptHs_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     ptHs_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

set xlabel 'y_{H_{soft}}'
plot yHs_ana u 5:6 w l t 'VBFNLO',\
     yHs_ten u (0.5*($1+$2)):(pbtofb*$3) w l t 'tensor',\
     yHs_ana u (0.5*($1+$2)):(pbtofb*$3) w l t 'analytic'

# set output 'comparison-vbfnlo-error.pdf'
# set xlabel 'p_{t,HH}'
# plot ptHH_ana u 5:(1.0):($7/$6) w e t 'VBFNLO',\
#      ptHH_ten u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'tensor',\
#      ptHH_ana u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'analytic'

# set xlabel 'y_{HH}'
# plot yHH_ana u 5:(1.0):($7/$6) w e t 'VBFNLO',\
#      yHH_ten u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'tensor',\
#      yHH_ana u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'analytic'

# set xlabel 'p_{t,H_1}'
# plot ptHa_ana u 5:(1.0):($7/$6) w e t 'VBFNLO',\
#      ptHa_ten u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'tensor',\
#      ptHa_ana u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'analytic'

# set xlabel 'y_{H_1}'
# plot yHa_ana u 5:(1.0):($7/$6) w e t 'VBFNLO',\
#      yHa_ten u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'tensor',\
#      yHa_ana u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'analytic'

# set xlabel 'p_{t,H_2}'
# plot ptHb_ana u 5:(1.0):($7/$6) w e t 'VBFNLO',\
#      ptHb_ten u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'tensor',\
#      ptHb_ana u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'analytic'

# set xlabel 'y_{H_2}'
# plot yHb_ana u 5:(1.0):($7/$6) w e t 'VBFNLO',\
#      yHb_ten u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'tensor',\
#      yHb_ana u (0.5*($1+$2)):(1.0):($4/$3) w e dt 2 t 'analytic'

set output 'comparison-vbfnlo-ratio.pdf'
set yrange [0.92:1.08]
set xlabel 'p_{t,HH}'
plot ptHH_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     ptHH_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'y_{HH}'
plot yHH_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     yHH_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'm_{HH}'
plot mHH_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     mHH_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel '{/Symbol D}{/Symbol f}_{HH}'
plot phiHH_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     phiHH_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'p_{t,H_1}'
plot ptHa_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     ptHa_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'y_{H_1}'
plot yHa_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     yHa_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'm_{H_1}'
plot mHa_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     mHa_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'p_{t,H_2}'
plot ptHb_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     ptHb_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'y_{H_2}'
plot yHb_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     yHb_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'm_{H_2}'
plot mHb_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     mHb_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'p_{t,H_{hard}}'
plot ptHd_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     ptHd_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'y_{H_{hard}}'
plot yHd_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     yHd_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'p_{t,H_{soft}}'
plot ptHs_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     ptHs_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

set xlabel 'y_{H_{soft}}'
plot yHs_ten u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'tensor/VBFNLO',\
     yHs_ana u (0.5*($1+$2)):(pbtofb*$3/$6) w l t 'analytic/VBFNLO'

