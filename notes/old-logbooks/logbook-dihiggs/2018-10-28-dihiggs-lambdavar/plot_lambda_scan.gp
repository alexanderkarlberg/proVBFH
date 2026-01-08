set term pdfcairo enhanced color transparent font "Times, 26" size 16cm,24cm
# set term wxt persist

fn = '<mergeidx.pl -f lambda_scan.dat '

set output 'lambda_scan.pdf'
set macros

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

pbtofb=1000
# plot cross section
reset
#set grid front
set logscale y
set xtics format ""
set ylabel '{/Symbol s} [fb]'
# set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
set label 2002 'HE-LHC 27 TeV' at graph 0.03,0.29
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.03,0.19
set label 2001 'Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.03,0.09
#set key maxrows 2
set key bottom right
set key width 1
set ytics 0.1,10,10000
set yrange [1:1000]
set xrange [-10:10]
set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'red' lw 3 dt 2
set multiplot
set origin 0.0,0.3
set size 1.0,0.7
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.57
plot fn.' LO.cross'   u 1:($3*pbtofb):($4*pbtofb) w @LOfill   t 'LO',\
     ''               u 1:($3*pbtofb)   w @LOlinethin not,\
     ''               u 1:($4*pbtofb)   w @LOlinethin not,\
     fn.' NLO.cross'  u 1:($3*pbtofb):($4*pbtofb) w @NLOfill  t 'NLO',\
     ''               u 1:($3*pbtofb)   w @NLOlinethin not,\
     ''               u 1:($4*pbtofb)   w @NLOlinethin not,\
     fn.' NNLO.cross' u 1:($3*pbtofb):($4*pbtofb) w @NNLOfill t 'NNLO',\
     ''               u 1:($3*pbtofb)   w @NNLOlinethin not,\
     ''               u 1:($4*pbtofb)   w @NNLOlinethin not,\
     fn.' N3LO.cross' u 1:($3*pbtofb):($4*pbtofb) w @N3LOfill t 'N^{3}LO',\
     fn.' LO.cross'   u 1:($2*pbtofb) w @LOline   not,\
     fn.' NLO.cross'  u 1:($2*pbtofb) w @NLOline  not,\
     fn.' NNLO.cross' u 1:($2*pbtofb) w @NNLOline not,\
     fn.' N3LO.cross' u 1:($2*pbtofb) w @N3LOline not

set yrange[0.97:1.07]
set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'red' lw 3 dt 2
set tmargin at screen 0.57
set bmargin at screen 0.15
set nologscale y
set title ''
set format x
set format y
set ytics 0.96,0.02,1.2
set ytics add ("" 1.08)
set mytics 2
unset label
set xlabel '{/Symbol k}' offset 0,0.5
set ylabel 'ratio to N^{3}LO'
plot fn.' LO.cross   N3LO.cross' u 1:($3/$6):($4/$6) w @LOfill   not,\
     ''               u 1:($3/$6)   w @LOlinethin not,\
     ''               u 1:($4/$6)   w @LOlinethin not,\
     fn.' NLO.cross  N3LO.cross' u 1:($3/$6):($4/$6) w @NLOfill  not,\
     ''               u 1:($3/$6)   w @NLOlinethin not,\
     ''               u 1:($4/$6)   w @NLOlinethin not,\
     fn.' NNLO.cross N3LO.cross' u 1:($3/$6):($4/$6) w @NNLOfill not,\
     ''               u 1:($3/$6)   w @NNLOlinethin not,\
     ''               u 1:($4/$6)   w @NNLOlinethin not,\
     fn.' N3LO.cross'            u 1:($3/$2):($4/$2) w @N3LOfill not,\
     fn.' LO.cross   N3LO.cross' u 1:($2/$6) w @LOline   not,\
     fn.' NLO.cross  N3LO.cross' u 1:($2/$6) w @NLOline  not,\
     fn.' NNLO.cross N3LO.cross' u 1:($2/$6) w @NNLOline not,\
     fn.' N3LO.cross'            u 1:($2/$2) w @N3LOline not

unset multiplot
set output
