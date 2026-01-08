set term pdfcairo enhanced color transparent font "Times, 26" size 20cm,12cm
# set term wxt persist

pretitle='{/*1.2 '
tinySpace='{/*0.5 &.}'
fn = 'n3lo_hist.dat'

colLO   = '#f4e0a6'
colNLO  = '#BFEFFF'
colNNLO = '#37BC61'
colN3LO = '#FF4500'

LOline   = "lines lc rgb '#f9B919' lw 2.0"
NLOline  = "lines lc rgb '#38B0DE' lw 2.0"
NNLOline = "lines lc rgb '#37BC61' lw 2.0"
NNLOlinethin = "lines lc rgb '#37BC61' lw 0.5"
N3LOline = "lines lc rgb '#FF4500' lw 2.0"

LOfill   = "filledcurves lc rgb '#FFF0C4' lw 2.0 fs transparent pattern 3"
NLOfill  = "filledcurves lc rgb '#05AEEFFF' lw 2.0 fs transparent solid 0.85"
NNLOfill = "filledcurves lc rgb '#37BC61' lw 2.0 fs transparent pattern 1"
N3LOfill = "filledcurves lc rgb '#FF4500' lw 2.0 fs transparent solid 0.5"

set output 'n3lo-total.pdf'

# plot ptH
reset
# set grid front
set xtics format ""
#set ytics 0.001,10,100
set yrange [3.85:4.2]
set ylabel '{/Symbol s} [pb]'
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.04,0.08
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.04,0.16
set label 2002 'LHC 13 TeV' at graph 0.04,0.23
# set label 1000 '{/*4 PRELIMINARY' at graph 0.45,0.50 center rotate by 25 tc rgb '#d0d0d0'

plot '<mergeidx.pl -f '.fn.' LO.cross'   u (3.) :(($2+$3)/2.):(2.5):(3.5):2:3 \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb colLO linewidth 2.0 ti 'LO',\
     '<mergeidx.pl -f '.fn.' NLO.cross'  u (3.) :(($2+$3)/2.):(4.5):(5.5):2:3 \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb colNLO linewidth 2.0 ti 'NLO',\
     '<mergeidx.pl -f '.fn.' NNLO.cross' u (3.) :(($2+$3)/2.):(6.5):(7.5):2:3 \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb colNNLO linewidth 2.0 ti 'NNLO',\
     '<mergeidx.pl -f '.fn.' N3LO.cross' u (3.) :(($2+$3)/2.):(8.5):(9.5):2:3 \
        w boxxyerrorbars fill transparent solid 0.5 lc rgb colN3LO linewidth 2.0 ti 'N3LO',\
     '<mergeidx.pl -f '.fn.' LO.cross'   u (3.) :(($2+$3)/2.):(2.5):(3.5) w xerrorbars ps 0 lc rgb colLO lw 3.0 not,\
     '<mergeidx.pl -f '.fn.' NLO.cross'  u (3.) :(($2+$3)/2.):(4.5):(5.5) w xerrorbars ps 0 lc rgb colNLO lw 3.0 not,\
     '<mergeidx.pl -f '.fn.' NNLO.cross' u (3.) :(($2+$3)/2.):(6.5):(7.5) w xerrorbars ps 0 lc rgb colNNLO lw 3.0 not,\
     '<mergeidx.pl -f '.fn.' N3LO.cross' u (3.) :(($2+$3)/2.):(8.5):(9.5) w xerrorbars ps 0 lc rgb colN3LO lw 3.0 not,\

reset

set output 'n3lo-hists.pdf'
set term pdfcairo size 15cm,22cm
# set grid front
# set label 1000 '{/*4 PRELIMINARY' at graph 0.45,0.50 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.04,0.08
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.04,0.16
set label 2002 'LHC 13 TeV' at graph 0.04,0.23
set key maxrows 2
set key top right
set key width 2
set key samplen 3
set format y "10^{%L}"
set xtics format ""
set ylabel ''
set xlabel ''
set title pretitle.'{d{/Symbol s}}/{d{p_{t,'.tinySpace.'H}}} [pb/GeV]' offset 0,-0.5

# set ylabel '{/Symbol s} [pb]'

set multiplot
set origin 0.0,0.3
set size 1.0,1.0
iindex = 7
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.47
set logscale y
set xrange [0:300]
plot '<mergeidx.pl -f '.fn.' LO.ptH' u (($1+$2)/2):4:5  w @LOfill ti 'LO',\
     '' u (($1+$2)/2):3 w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH' u (($1+$2)/2.):4:5  w @NLOfill ti 'NLO',\
     '' u (($1+$2)/2):3 w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH' u (($1+$2)/2.):4:5 w @NNLOfill ti 'NNLO',\
     '' u (($1+$2)/2):4 w @NNLOlinethin not,\
     '' u (($1+$2)/2):5 w @NNLOlinethin not,\
     '' u (($1+$2)/2):3 w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH' u (($1+$2)/2.):4:5 w @N3LOfill ti 'N3LO',\
     '' u (($1+$2)/2):4  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):5 w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):3 w @N3LOline not,\


set tmargin at screen 0.47
set bmargin at screen 0.115
set nologscale y
set title ''
set ytics 0.9,0.01,1.025
set mytics 2
set format x
set format y
unset label
#set ylabel 'ratio to N3LO'
set xlabel 'p_{t,H} [GeV]'
set yrange [0.99:1.03]
set xrange [0:300]
plot '<mergeidx.pl -f '.fn.' LO.ptH N3LO.ptH' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH N3LO.ptH' u (($1+$2)/2.):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH N3LO.ptH' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '' u (($1+$2)/2):($4/$8) w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8) w @NNLOlinethin not,\
     '' u (($1+$2)/2.):($4/$8):($5/$8) w @NNLOfill not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH' u (($1+$2)/2):($3/$3) w @N3LOline not,\
     '' u (($1+$2)/2):($4/$3)  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):($5/$3) w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\

unset multiplot

# YH plot
# set grid
# set label 1000 '{/*4 PRELIMINARY' at graph 0.45,0.50 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.04,0.08
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.04,0.16
set label 2002 'LHC 13 TeV' at graph 0.04,0.23
set key maxrows 2
set key top right
set key width 2
set key samplen 3
# set format y "10^{%L}"
set xtics format ""
set xlabel ''
set ylabel ''
set title pretitle.'{d{/Symbol s}}/{d{y_{H}}} [pb]' offset 0,-0.5
# set title 'd{/Symbol s}/d y_H [pb]'
set multiplot
set origin 0.0,0.3
set size 1.0,1.0
iindex = 7
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.47
set nologscale y
set mytics default
set ytics 0,0.2,1.5
set mytics 2
# set ytics 0.0001,10,100
set yrange [0:1.2]
set xrange [0:4.25]
plot '<mergeidx.pl -f '.fn.' LO.yH'   u (($1+$2)/2):3 w @LOline not,\
     '' u (($1+$2)/2):4:5  w @LOfill ti 'LO',\
     '<mergeidx.pl -f '.fn.' NLO.yH'  u (($1+$2)/2):3 w @NLOline not,\
     '' u (($1+$2)/2.):4:5  w @NLOfill ti 'NLO',\
     '<mergeidx.pl -f '.fn.' NNLO.yH' u (($1+$2)/2):3 w @NNLOline not,\
     '' u (($1+$2)/2):4 w @NNLOlinethin not,\
     '' u (($1+$2)/2):5 w @NNLOlinethin not,\
     '' u (($1+$2)/2.):4:5 w @NNLOfill ti 'NNLO',\
     '<mergeidx.pl -f '.fn.' N3LO.yH' u (($1+$2)/2):3 w @N3LOline not,\
     '' u (($1+$2)/2):4  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):5 w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2.):4:5 w @N3LOfill ti 'N3LO'

set tmargin at screen 0.47
set bmargin at screen 0.115
set nologscale y
set title ''
set ytics 0.9,0.01,1.025
set mytics 2
set format x
set format y
#set ylabel 'ratio to N3LO'
unset label
set xlabel 'y_{H}'
set yrange [0.99:1.03]
set xrange [0:4.25]
plot '<mergeidx.pl -f '.fn.' LO.yH N3LO.yH' u (($1+$2)/2):($4/$8):($5/$8)  w @LOfill not,\
     '' u (($1+$2)/2):($3/$8) w @LOline not,\
     '<mergeidx.pl -f '.fn.' NLO.yH N3LO.yH' u (($1+$2)/2.):($4/$8):($5/$8)  w @NLOfill not,\
     '' u (($1+$2)/2):($3/$8) w @NLOline not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH N3LO.yH' u (($1+$2)/2.):($4/$8):($5/$8) w @NNLOfill not,\
     '' u (($1+$2)/2):($4/$8) w @NNLOlinethin not,\
     '' u (($1+$2)/2):($5/$8) w @NNLOlinethin not,\
     '' u (($1+$2)/2):($3/$8) w @NNLOline not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH'  u (($1+$2)/2.):($4/$3):($5/$3) w @N3LOfill not,\
     '' u (($1+$2)/2):($4/$3)  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):($5/$3) w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):($3/$3) w @N3LOline not,\

unset multiplot

set output
