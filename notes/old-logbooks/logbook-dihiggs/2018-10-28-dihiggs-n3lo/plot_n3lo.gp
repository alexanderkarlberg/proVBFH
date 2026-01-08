set term pdfcairo enhanced color transparent font "Times, 26" size 14cm,24cm
# set term wxt persist

fn = 'n3lo_hist_27tev.dat'
rbn=' | rebin.pl -1 -r -c 5'


LOline   = "lines lc rgb '#f9B919' lw 2.0"
LOlinethin   = "lines lc rgb '#fccc55' lw 0.5"
NLOline  = "lines lc rgb '#38B0DE' lw 2.0"
NLOlinethin  = "lines lc rgb '#05AEEFFF' lw 0.5"
NNLOline = "lines lc rgb '#37BC61' lw 2.0"
NNLOlinethin = "lines lc rgb '#37BC61' lw 0.5"
N3LOline = "lines lc rgb '#FF4500' lw 2.0"

LOfill   = "filledcurves lc rgb '#fccc55' lw 2.0 fs transparent pattern 4"
#LOfill   = "filledcurves lc rgb '#FFF0C4' lw 2.0 fs transparent solid 0.8"
NLOfill  = "filledcurves lc rgb '#05AEEFFF' lw 2.0 fs transparent pattern 5"
# NLOfill  = "filledcurves lc rgb '#05AEEFFF' lw 2.0 fs transparent solid 0.85"
NNLOfill = "filledcurves lc rgb '#37BC61' lw 2.0 fs transparent pattern 1"
N3LOfill = "filledcurves lc rgb '#FF4500' lw 2.0 fs transparent solid 0.5"

set output 'n3lo_27tev.pdf'

pbtofb=1000
set title ''

# plot ptHH
reset
set xtics format ""
set ylabel '{/Symbol s} [fb]'
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.97,0.04 right
set label 2001 'Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.10 right
set label 2002 'HE-LHC 27 TeV' at graph 0.97,0.16 right

plot '<mergeidx.pl -f '.fn.' LO.cross'   u (3.) :(pbtofb*($2+$3)/2.):(2.5):(3.5):(pbtofb*$2):(pbtofb*$3) \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb '#FFF8DC' linewidth 2.0 ti 'LO',\
     '<mergeidx.pl -f '.fn.' NLO.cross'  u (3.) :(pbtofb*($2+$3)/2.):(4.5):(5.5):(pbtofb*$2):(pbtofb*$3) \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb '#BFEFFF' linewidth 2.0 ti 'NLO',\
     '<mergeidx.pl -f '.fn.' NNLO.cross' u (3.) :(pbtofb*($2+$3)/2.):(6.5):(7.5):(pbtofb*$2):(pbtofb*$3) \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb '#37BC61' linewidth 2.0 ti 'NNLO',\
     '<mergeidx.pl -f '.fn.' N3LO.cross' u (3.) :(pbtofb*($2+$3)/2.):(8.5):(9.5):(pbtofb*$2):(pbtofb*$3) \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb '#FF4500' linewidth 2.0 ti 'N3LO',\
     '<mergeidx.pl -f '.fn.' LO.cross'   u (3.) :(pbtofb*($2+$3)/2.):(2.5):(3.5) w xerrorbars ps 0 lc rgb '#FFF8DC' lw 3.0 not,\
     '<mergeidx.pl -f '.fn.' NLO.cross'  u (3.) :(pbtofb*($2+$3)/2.):(4.5):(5.5) w xerrorbars ps 0 lc rgb '#BFEFFF' lw 3.0 not,\
     '<mergeidx.pl -f '.fn.' NNLO.cross' u (3.) :(pbtofb*($2+$3)/2.):(6.5):(7.5) w xerrorbars ps 0 lc rgb '#37BC61' lw 3.0 not,\
     '<mergeidx.pl -f '.fn.' N3LO.cross' u (3.) :(pbtofb*($2+$3)/2.):(8.5):(9.5) w xerrorbars ps 0 lc rgb '#FF4500' lw 3.0 not,\

reset
#set grid
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dp_{t,HH} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.94 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.86 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.03,0.94
set yrange [0.003:0.1]
set xtics format ""
set logscale y
set key maxrow 2 width 2
set key bottom left
set ytics 0.001,10,100
plot '<mergeidx.pl -f '.fn.' LO.ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptHH' u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.ptHH'   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptHH'  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptHH' u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptHH' u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115
set yrange [0.97:1.03]
set xrange [0:300]
set nologscale y
set format x
unset label 2000
unset label 2001
unset label 2002
set ytics 0.9,0.01,1.02
set ylabel 'ratio to N^{3}LO' offset 1,0
set xlabel 'p_{t,HH} [GeV]' offset 0,0.3
plot '<mergeidx.pl -f '.fn.' LO.ptHH N3LO.ptHH' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.ptHH N3LO.ptHH' u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptHH N3LO.ptHH' u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptHH' u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.ptHH N3LO.ptHH'   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptHH N3LO.ptHH'  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptHH N3LO.ptHH' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptHH' u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set xlabel ''
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dp_{t,H_1} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.94 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.86 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.03,0.94
set yrange [0.002:0.1]
set xtics format ""
set logscale y
set ytics 0.001,10,100
plot '<mergeidx.pl -f '.fn.' LO.ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_h' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_h' u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.ptH_h'   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_h'  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_h' u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_h' u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115
set xlabel 'p_{t,H_1} [GeV]' offset 0,0.3
set ylabel 'ratio to N^{3}LO' offset 1,0
set ytics 0.9,0.01,1.02
set nologscale y
set format x
set yrange [0.97:1.03]
unset label 2000
unset label 2001
unset label 2002
set xrange [0:300]
plot '<mergeidx.pl -f '.fn.' LO.ptH_h N3LO.ptH_h' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_h N3LO.ptH_h' u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_h N3LO.ptH_h' u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_h' u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.ptH_h N3LO.ptH_h'   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_h N3LO.ptH_h'  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_h N3LO.ptH_h' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_h' u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set xlabel ''
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dp_{t,H_2} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.94 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.86 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.03,0.94
set yrange [0.005:0.12]
set xrange [0:200]
set xtics format ""
set logscale y
set ytics 0.001,10,100
plot '<mergeidx.pl -f '.fn.' LO.ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_s' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_s' u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.ptH_s'   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_s'  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_s' u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_s' u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115
set xlabel 'p_{t,H_2} [GeV]' offset 0,0.3
set ytics 0.9,0.01,1.02
set nologscale y
set ylabel 'ratio to N^{3}LO' offset 1,0
set format x
unset label 2000
unset label 2001
unset label 2002
set yrange [0.97:1.03]
set xrange [0:200]
plot '<mergeidx.pl -f '.fn.' LO.ptH_s N3LO.ptH_s' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_s N3LO.ptH_s' u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_s N3LO.ptH_s' u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_s' u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.ptH_s N3LO.ptH_s'   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH_s N3LO.ptH_s'  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH_s N3LO.ptH_s' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH_s' u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set xlabel ''
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dy_{HH} [fb]' #offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.86 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.78 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.97,0.94 right
set yrange [0.0:3.7]
set key maxrow 4
set key top left Left invert reverse
set xrange [-4.2:4.2]
set xtics format ""
#set logscale y
set ytics 0.0,0.5,10
plot '<mergeidx.pl -f '.fn.' LO.yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.yHH'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.yHH'.rbn u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.yHH'.rbn   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yHH'.rbn  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yHH'.rbn u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yHH'.rbn u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115

set yrange [0.97:1.03]
set ytics 0.9,0.01,1.02
set nologscale y
set format x
set ylabel 'ratio to N^{3}LO' offset 1,0
unset label 2000
unset label 2001
unset label 2002
set xlabel 'y_{HH}' offset 0,0.3
plot '<mergeidx.pl -f '.fn.' LO.yHH N3LO.yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.yHH N3LO.yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.yHH N3LO.yHH'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.yHH'.rbn u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.yHH N3LO.yHH'.rbn   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yHH N3LO.yHH'.rbn  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yHH N3LO.yHH'.rbn u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yHH'.rbn u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set xlabel ''
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dy_{H_1} [fb]' #offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.86 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.78 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.97,0.94 right
set yrange [0.0:2.6]
set xrange [-4.2:4.2]
set xtics format ""
#set logscale y
set ytics 0,0.5,10
plot '<mergeidx.pl -f '.fn.' LO.yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_h'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_h'.rbn u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.yH_h'.rbn   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_h'.rbn  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_h'.rbn u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_h'.rbn u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set yrange [0.97:1.03]
set ytics 0.9,0.01,1.02
set nologscale y
set format x
set tmargin at screen 0.55
set bmargin at screen 0.115
set xlabel 'y_{H_1}' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to N^{3}LO' offset 1,0
plot '<mergeidx.pl -f '.fn.' LO.yH_h N3LO.yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_h N3LO.yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_h N3LO.yH_h'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_h'.rbn u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.yH_h N3LO.yH_h'.rbn   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_h N3LO.yH_h'.rbn  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_h N3LO.yH_h'.rbn u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_h'.rbn u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set xlabel ''
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dy_{H_2} [fb]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.86 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.78 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.97,0.94 right
set yrange [0:2.5]
set xrange [-4.2:4.2]
set xtics format ""
#set logscale y
set ytics 0,0.5,10
plot '<mergeidx.pl -f '.fn.' LO.yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_s'.rbn u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_s'.rbn u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.yH_s'.rbn   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_s'.rbn  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_s'.rbn u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_s'.rbn u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115


set yrange [0.97:1.03]
set ytics 0.9,0.01,1.02
set nologscale y
set format x
set xlabel 'y_{H_2}' offset 0,0.3
unset label 2000
unset label 2001
unset label 2002
set ylabel 'ratio to N^{3}LO' offset 1,0
plot '<mergeidx.pl -f '.fn.' LO.yH_s N3LO.yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_s N3LO.yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_s N3LO.yH_s'.rbn u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_s'.rbn u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.yH_s N3LO.yH_s'.rbn   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yH_s N3LO.yH_s'.rbn  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH_s N3LO.yH_s'.rbn u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH_s'.rbn u (($1+$2)/2):($3/$3) w @N3LOline not

unset multiplot



reset
#set grid
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/dm_{HH} [fb/GeV]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.94 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.86 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.03,0.94
set yrange [0.0002:0.03]
set xrange [300:1800]
set xtics format ""
set logscale y
set key maxrow 2 width 2
set key bottom left
set ytics 0.0001,10,100
set xtics 0,300,2000
set mxtics 3
plot '<mergeidx.pl -f '.fn.' LO.mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.mHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.mHH' u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.mHH'   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.mHH'  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.mHH' u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.mHH' u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115
set yrange [0.97:1.03]
set xrange [300:1800]
set nologscale y
set format x
unset label 2000
unset label 2001
unset label 2002
set ytics 0.9,0.01,1.02
set ylabel 'ratio to N^{3}LO' offset 1,0
set xlabel 'm_{HH} [GeV]' offset 0,0.3
plot '<mergeidx.pl -f '.fn.' LO.mHH N3LO.mHH' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.mHH N3LO.mHH' u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.mHH N3LO.mHH' u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.mHH' u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.mHH N3LO.mHH'   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.mHH N3LO.mHH'  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.mHH N3LO.mHH' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.mHH' u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set xlabel ''
set multiplot
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.55
set format x
set ylabel 'd{/Symbol s}/d{/Symbol f}_{HH} [fb]' offset 1,0
set label 2000 '{/*0.8PDF4LHC15\_nnlo\_mc' at graph 0.97,0.08 right
set label 2001 '{/*0.8Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.97,0.17 right
set label 2002 '{/*0.8HE-LHC 27 TeV' at graph 0.97,0.26 right
set xrange [0:3.15]
set xtics format ""
set logscale y
set yrange [1:8]
set key top left
set key maxrow 4
set ytics 0.001,10,1000
set ytics add (2,3,5,8,10)
plot '<mergeidx.pl -f '.fn.' LO.phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @LOfill t 'LO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NLOfill t 'NLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.phiHH' u (($1+$2)/2):($4*pbtofb):($5*pbtofb)  w @NNLOfill t 'NNLO',\
     '' u (($1+$2)/2):($4*pbtofb)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5*pbtofb)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.phiHH' u (($1+$2)/2.):($4*pbtofb):($5*pbtofb) w @N3LOfill t 'N^{3}LO',\
     '<mergeidx.pl -f '.fn.' LO.phiHH'   u (($1+$2)/2):($3*pbtofb) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.phiHH'  u (($1+$2)/2):($3*pbtofb) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.phiHH' u (($1+$2)/2):($3*pbtofb) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.phiHH' u (($1+$2)/2):($3*pbtofb) w @N3LOline not

set tmargin at screen 0.55
set bmargin at screen 0.115

set yrange [0.97:1.03]
set ytics 0.9,0.01,1.02
set nologscale y
set format x
set ylabel 'ratio to N^{3}LO' offset 1,0
unset label 2000
unset label 2001
unset label 2002
set xlabel '{/Symbol f}_{HH}' offset 0,0.3
plot '<mergeidx.pl -f '.fn.' LO.phiHH N3LO.phiHH' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @LOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @LOlinethin not,\
     '<mergeidx.pl -f '.fn.' NLO.phiHH N3LO.phiHH' u (($1+$2)/2):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NLOlinethin not,\
     '<mergeidx.pl -f '.fn.' NNLO.phiHH N3LO.phiHH' u (($1+$2)/2):($4/$8):($5/$8)  w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8)  w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8)  w @NNLOlinethin not,\
     '<mergeidx.pl -f '.fn.' N3LO.phiHH' u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '<mergeidx.pl -f '.fn.' LO.phiHH N3LO.phiHH'   u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.phiHH N3LO.phiHH'  u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.phiHH N3LO.phiHH' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.phiHH' u (($1+$2)/2):($3/$3) w @N3LOline not
unset multiplot

set output
