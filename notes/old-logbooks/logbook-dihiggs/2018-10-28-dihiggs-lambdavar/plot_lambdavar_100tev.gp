set term pdfcairo enhanced color transparent font "Times, 26" size 14cm,22cm
#set term pdfcairo enhanced color transparent font "Times, 26" size 13cm,18cm
# set term wxt persist
fn100 = '-f 100tev_n3lo_hist_lambda100.dat '
fn080 = '-f 100tev_n3lo_hist_lambda080.dat '
fn090 = '-f 100tev_n3lo_hist_lambda090.dat '
fn110 = '-f 100tev_n3lo_hist_lambda110.dat '
fn120 = '-f 100tev_n3lo_hist_lambda120.dat '
fn100_LO = '-f 100tev_lo_hist_lambda100.dat '
fn080_LO = '-f 100tev_lo_hist_lambda080.dat '
fn090_LO = '-f 100tev_lo_hist_lambda090.dat '
fn110_LO = '-f 100tev_lo_hist_lambda110.dat '
fn120_LO = '-f 100tev_lo_hist_lambda120.dat '
rbn=' | rebin.pl -1 -r -c 5'
pbtofb=1000

set output 'lambdavar_100tev.pdf'

# plot ptHH
reset
#set grid
# set xtics format ""
# set ylabel '{/Symbol s}'

# plot '<mergeidx.pl -f '.fn.' LO.cross'   u (3.) :(($2+$3)/2.):(2.5):(3.5):2:3 \
#         w boxxyerrorbars fill transparent solid 0.5 lc rgb '#FFF8DC' linewidth 2.0 ti 'LO',\
#      '<mergeidx.pl -f '.fn.' NLO.cross'  u (3.) :(($2+$3)/2.):(4.5):(5.5):2:3 \
#         w boxxyerrorbars fill transparent solid 0.5 lc rgb '#BFEFFF' linewidth 2.0 ti 'NLO',\
#      '<mergeidx.pl -f '.fn.' NNLO.cross' u (3.) :(($2+$3)/2.):(6.5):(7.5):2:3 \
#         w boxxyerrorbars fill transparent solid 0.5 lc rgb '#37BC61' linewidth 2.0 ti 'NNLO',\
#      '<mergeidx.pl -f '.fn.' N3LO.cross' u (3.) :(($2+$3)/2.):(8.5):(9.5):2:3 \
#         w boxxyerrorbars fill transparent solid 0.5 lc rgb '#FF4500' linewidth 2.0 ti 'N3LO',\
#      '<mergeidx.pl -f '.fn.' LO.cross'   u (3.) :(($2+$3)/2.):(2.5):(3.5) w xerrorbars ps 0 lc rgb '#FFF8DC' lw 3.0 not,\
#      '<mergeidx.pl -f '.fn.' NLO.cross'  u (3.) :(($2+$3)/2.):(4.5):(5.5) w xerrorbars ps 0 lc rgb '#BFEFFF' lw 3.0 not,\
#      '<mergeidx.pl -f '.fn.' NNLO.cross' u (3.) :(($2+$3)/2.):(6.5):(7.5) w xerrorbars ps 0 lc rgb '#37BC61' lw 3.0 not,\
#      '<mergeidx.pl -f '.fn.' N3LO.cross' u (3.) :(($2+$3)/2.):(8.5):(9.5) w xerrorbars ps 0 lc rgb '#FF4500' lw 3.0 not,\

# reset
# set grid
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dp_{t,HH} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.03:1.1]
set xtics format ""
set logscale y
set ytics 0.001,10,100
set xrange [0:300]
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
plot '<mergeidx.pl '.fn080.fn100.'ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'ptHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'ptHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'ptHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'ptHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'ptHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set xrange [0:300]
set yrange [0.8:1.25]
set nologscale y
set xlabel 'p_{t,HH} [GeV]' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set format x
set ytics 0.8,0.1,1.2
set key bottom right reverse samplen 1.3 width 0.5
#set key at 1.5,0.0123
plot '<mergeidx.pl '.fn080.fn100.'ptHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'ptHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'ptHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'ptHH' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'ptHH' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'ptHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'

unset multiplot

set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dp_{t,H_1} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.03:1]
set xtics format ""
set logscale y
set ytics 0.001,10,100
set xrange [0:300]
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
plot '<mergeidx.pl '.fn080.fn100.'ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'ptH_h'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'ptH_h'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'ptH_h'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'ptH_h'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'ptH_h'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set xrange [0:300]
set yrange [0.8:1.25]
set nologscale y
set format x
set xlabel 'p_{t,H_1} [GeV]' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set ytics 0.8,0.1,1.2
set key bottom left reverse samplen 1.3 width 0.5
#set key at 1.5,0.0123
plot '<mergeidx.pl '.fn080.fn100.'ptH_h' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'ptH_h' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'ptH_h' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'ptH_h' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'ptH_h' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'ptH_h'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'
unset multiplot


set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dp_{t,H_2} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.01:1.5]
set xtics format ""
set logscale y
set ytics 0.001,10,100
set xrange [0:300]
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
plot '<mergeidx.pl '.fn080.fn100.'ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'ptH_s'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'ptH_s'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'ptH_s'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'ptH_s'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'ptH_s'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set xrange [0:300]
set yrange [0.8:1.25]
set nologscale y
set format x
set xlabel 'p_{t,H_2} [GeV]' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set ytics 0.8,0.1,1.2
set key bottom left reverse samplen 1.3 width 0.5
#set key at 1.5,0.0123
plot '<mergeidx.pl '.fn080.fn100.'ptH_s' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'ptH_s' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'ptH_s' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'ptH_s' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'ptH_s' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'ptH_s'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'
unset multiplot

set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dy_{HH} [fb]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.1:100]
set xtics format ""
set logscale y
set ytics 0.001,10,100
set xrange [-4.2:4.2]
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
set key at -1.5,0.123
plot '<mergeidx.pl '.fn080.fn100.'yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'yHH'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'yHH'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'yHH'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'yHH'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'yHH'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set yrange [0.8:1.25]
set nologscale y
set format x
set ytics 0.8,0.1,1.2
set xrange [-4.2:4.2]
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set key bottom left reverse samplen 1.3 width 0.5
set key at 2.2,1.007
set xlabel 'y_{HH}' offset 0,0.3
plot '<mergeidx.pl '.fn080.fn100.'yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'yHH'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'

unset multiplot



set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dy_{H_1} [fb]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.1:100]
set xtics format ""
set logscale y
set ytics 0.001,10,100
set xrange [-4.2:4.2]
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
set key at -1.5,0.123
plot '<mergeidx.pl '.fn080.fn100.'yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'yH_h'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'yH_h'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'yH_h'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'yH_h'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'yH_h'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set yrange [0.8:1.25]
set nologscale y
set format x
set ytics 0.8,0.1,1.2
set xrange [-4.2:4.2]
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set xlabel 'y_{H_1}' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set key bottom left reverse samplen 1.5 width -1.2
set key at -1.5,0.123
plot '<mergeidx.pl '.fn080.fn100.'yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'yH_h'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'


set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dy_{H_2} [fb]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.1:100]
set xtics format ""
set logscale y
set ytics 0.001,10,100
set xrange [-4.2:4.2]
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
set key at -1.5,0.123
plot '<mergeidx.pl '.fn080.fn100.'yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'yH_s'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'yH_s'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'yH_s'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'yH_s'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'yH_s'.rbn   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set yrange [0.8:1.25]
set nologscale y
set format x
set ytics 0.8,0.1,1.2
set xrange [-4.2:4.2]
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set xlabel 'y_{H_2}' offset 0,0.3
set key bottom left reverse samplen 1.3 width 0.5
set key at -3.79,0.8118
unset label 2000
unset label 2001
unset label 2002
plot '<mergeidx.pl '.fn080.fn100.'yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'yH_s'.rbn   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'
unset multiplot

reset
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dm_{HH} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.93 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.84 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.02,0.93
set yrange [0.002:0.3]
set xtics format ""
set logscale y
set ytics 0.0001,10,100
set xtics 0,300,2000
set xrange [300:1800]
set mxtics 3
set xlabel ''
set key bottom left reverse samplen 1.5 width -1.2
plot '<mergeidx.pl '.fn080.fn100.'mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'mHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'mHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'mHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'mHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'mHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set xrange [300:1800]
set yrange [0.8:1.25]
set nologscale y
set xlabel 'm_{HH} [GeV]' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set format x
set key bottom right reverse samplen 1.3 width 0.5
set key at 1760,0.8118
set ytics 0.8,0.1,1.2
plot '<mergeidx.pl '.fn080.fn100.'mHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'mHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'mHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'mHH' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'mHH' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'mHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'
unset multiplot

reset
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/d{/Symbol f}_{HH} [fb]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.98,0.07 right
set label 2001 '{/*0.8Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.98,0.16 right
set label 2002 '{/*0.8FCC 100 TeV' at graph 0.98,0.25 right
set yrange [10:100]
set xtics format ""
set logscale y
set ytics 0.001,10,1000
set ytics add (2,3,5,8,10)
set xrange [0:3.15]
set xlabel ''
set key top left reverse samplen 1.5 width -1.2
plot '<mergeidx.pl '.fn080.fn100.'phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.8',\
     '<mergeidx.pl '.fn090.fn100.'phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 0.9',\
     '<mergeidx.pl '.fn100.fn100.'phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1   ',\
     '<mergeidx.pl '.fn110.fn100.'phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.1',\
     '<mergeidx.pl '.fn120.fn100.'phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 t '{/Symbol k} = 1.2',\
     '<mergeidx.pl '.fn080.fn100.'phiHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'phiHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'phiHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'phiHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'phiHH'   u (($1+$2)/2):($3*pbtofb) w lines lw 2.0 lc 3 not

set tmargin at screen 0.55
set bmargin at screen 0.115
set yrange [0.8:1.25]
set nologscale y
set format x
set ytics 0.8,0.1,1.2
set xrange [0:3.15]
unset label 2000
unset label 2001
unset label 2002
set key bottom left reverse samplen 1.3 width 0.5
set key at 0.1,0.8118
set ylabel 'ratio to {/Symbol k}=1' offset 1,0
set xlabel '{/Symbol f}_{HH}' offset 0,0.3
plot '<mergeidx.pl '.fn080.fn100.'phiHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 1 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn090.fn100.'phiHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 2 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn100.fn100.'phiHH' u (($1+$2)/2):($4/$8):($5/$8) \
     w filledcurves lc 7 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn110.fn100.'phiHH' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 4 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn120.fn100.'phiHH' u (($1+$2)/2):($4/$8):($5/$8)\
     w filledcurves lc 3 lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl '.fn080_LO.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 1 not,\
     '<mergeidx.pl '.fn090_LO.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 2 not,\
     '<mergeidx.pl '.fn100_LO.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 7 not,\
     '<mergeidx.pl '.fn110_LO.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 4 not,\
     '<mergeidx.pl '.fn120_LO.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 dt 2 lc 3 not,\
     '<mergeidx.pl '.fn080.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 1 not,\
     '<mergeidx.pl '.fn090.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 2 not,\
     '<mergeidx.pl '.fn100.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 7 not,\
     '<mergeidx.pl '.fn110.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 4 not,\
     '<mergeidx.pl '.fn120.fn100.'phiHH'   u (($1+$2)/2):($3/$8) w lines lw 2.0 lc 3 not,\
     sqrt(-1) w l lw 2.0 lc rgb 'black' dt 2 t 'LO'

unset multiplot


set output
