set term pdfcairo enhanced color transparent font "Times, 26" size 22cm,20cm
# set term wxt persist

fn = '<mergeidx.pl -f n3lo_energy-scan.dat '

set output 'n3lo_energy.pdf'
set macros

LOline   = "lines lc rgb '#f9B919' lw 2.0"
NLOline  = "lines lc rgb '#38B0DE' lw 2.0"
NNLOline = "lines lc rgb '#37BC61' lw 2.0"
NNLOlinethin = "lines lc rgb '#37BC61' lw 0.5"
N3LOline = "lines lc rgb '#FF4500' lw 2.0"

LOfill   = "filledcurves lc rgb '#FFF0C4' lw 2.0 fs transparent pattern 3"
NLOfill  = "filledcurves lc rgb '#05AEEFFF' lw 2.0 fs transparent solid 0.85"
NNLOfill = "filledcurves lc rgb '#37BC61' lw 2.0 fs transparent pattern 1"
N3LOfill = "filledcurves lc rgb '#FF4500' lw 2.0 fs transparent solid 0.5"

# plot cross section
reset
set grid front
set logscale xy
set xtics format ""
set ylabel '{/Symbol s} [pb]'
# set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'NNPDF30\_nnlo\_as\_0118' at graph 0.02,0.92
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.02,0.83
set key maxrows 2
set key bottom right
set key width 3
set ytics 0.1,10,100
set ytics add ("" 0.01)
set xtics add (7,13,20,30)
set yrange [1:100]
set xrange [7:100]
set multiplot
set origin 0.0,0.3
set size 1.0,0.7
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.57
plot fn.' LO.cross'   u 1:3:4 w @LOfill   t 'LO',\
     fn.' NLO.cross'  u 1:3:4 w @NLOfill  t 'NLO',\
     fn.' NNLO.cross' u 1:3:4 w @NNLOfill t 'NNLO',\
     ''               u 1:3   w @NNLOlinethin not,\
     ''               u 1:4   w @NNLOlinethin not,\
     fn.' N3LO.cross' u 1:3:4 w @N3LOfill t 'N3LO',\
     fn.' LO.cross'   u 1:2 w @LOline   not,\
     fn.' NLO.cross'  u 1:2 w @NLOline  not,\
     fn.' NNLO.cross' u 1:2 w @NNLOline not,\
     fn.' N3LO.cross' u 1:2 w @N3LOline not

set yrange[0.98:1.08]
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
set xlabel '{/Symbol \326} s [TeV]' offset 0,0.5
set ylabel 'ratio to N3LO'
plot fn.' LO.cross   N3LO.cross' u 1:($3/$6):($4/$6) w @LOfill   not,\
     fn.' NLO.cross  N3LO.cross' u 1:($3/$6):($4/$6) w @NLOfill  not,\
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
