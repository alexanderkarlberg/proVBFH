set term pdfcairo enhanced color transparent font "Times, 26" size 20cm,15cm
# set term wxt persist

fn = '<mergeidx.pl -f n3lo_energy-scan.dat '

set output 'n3lo_energy.pdf'
set macros

LOline   = "lines lc rgb '#FFCC11' lw 2.0"
NLOline  = "lines lc rgb '#0099CC' lw 2.0"
NNLOline = "lines lc rgb '#37BC61' lw 2.0"
N3LOline = "lines lc rgb '#FF4500' lw 2.0"

LOfill   = "filledcurves lc rgb '#FFCC11' lw 2.0 fs transparent solid 0.3"
NLOfill  = "filledcurves lc rgb '#0099CC' lw 2.0 fs transparent solid 0.3"
NNLOfill = "filledcurves lc rgb '#37BC61' lw 2.0 fs transparent solid 0.3"
N3LOfill = "filledcurves lc rgb '#FF4500' lw 2.0 fs transparent solid 0.3"

# plot cross section
reset
set grid
set ylabel '{/Symbol s} [pb]'
set xlabel '{/Symbol \326} s [TeV]'
set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'NNPDF30\_nnlo\_as\_118' at graph 0.02,0.94
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.02,0.88
set key maxrows 2
set key bottom right
set key width 3
plot fn.' LO.cross'   u 1:3:4 w @LOfill   t 'LO',\
     fn.' NLO.cross'  u 1:3:4 w @NLOfill  t 'NLO',\
     fn.' NNLO.cross' u 1:3:4 w @NNLOfill t 'NNLO',\
     fn.' N3LO.cross' u 1:3:4 w @N3LOfill t 'N3LO',\
     fn.' LO.cross'   u 1:2 w @LOline   not,\
     fn.' NLO.cross'  u 1:2 w @NLOline  not,\
     fn.' NNLO.cross' u 1:2 w @NNLOline not,\
     fn.' N3LO.cross' u 1:2 w @N3LOline not

set yrange[0.95:1.1]
set ylabel 'ratio to N3LO'
plot fn.' LO.cross   N3LO.cross' u 1:($3/$6):($4/$6) w @LOfill t 'LO',\
     fn.' NLO.cross  N3LO.cross' u 1:($3/$6):($4/$6) w @NLOfill t 'NLO',\
     fn.' NNLO.cross N3LO.cross' u 1:($3/$6):($4/$6) w @NNLOfill t 'NNLO',\
     fn.' N3LO.cross'            u 1:($3/$2):($4/$2) w @N3LOfill t 'N3LO',\
     fn.' LO.cross   N3LO.cross' u 1:($2/$6) w @LOline   not,\
     fn.' NLO.cross  N3LO.cross' u 1:($2/$6) w @NLOline  not,\
     fn.' NNLO.cross N3LO.cross' u 1:($2/$6) w @NNLOline not,\
     fn.' N3LO.cross'            u 1:($2/$2) w @N3LOline not

set output
