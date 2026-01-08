set term pdfcairo enhanced color transparent font "Times, 26" size 16cm,15cm
# set term wxt persist

fn = 'n3lo_scale-scan.dat'

set output 'n3lo_scale-scan_thin.pdf'
set macros

LOline   = "lines lc rgb '#FFCC11' lw 2.0"
NLOline  = "lines lc rgb '#0099CC' lw 2.0"
NNLOline = "lines lc rgb '#37BC61' lw 2.0"
N3LOline = "lines lc rgb '#FF4500' lw 2.0"

# plot cross section
reset
set grid
set ylabel '{/Symbol s} [pb]'
set xlabel '{/Symbol x}'
#set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'NNPDF30\_nnlo\_as\_118' at graph 0.98,0.81 right
set label 2001 '{/Symbol m}={/Symbol m}_R={/Symbol m}_F={/Symbol x} m_h' at graph 0.98,0.88 right
set label 2002 'LHC 13 TeV' at graph 0.98,0.94 right
set key maxrows 2
set key bottom left Left reverse
set key width 3
set xrange [0.125:8]
# set xrange [0.1:10]
set xtics add (1/8.,1/4.,1/2.,2,4,8)
set log x
set yrange [3.8:4.15]
plot fn u 1:2 w @LOline t 'LO',\
     fn u 1:3 w @NLOline t 'NLO',\
     fn u 1:4 w @NNLOline t 'NNLO',\
     fn u 1:5 w @N3LOline t 'N3LO'

set yrange [0.95:1.05]
set ylabel 'ratio to N3LO'
plot fn u 1:($2/$5) w @LOline t 'LO',\
     fn u 1:($3/$5) w @NLOline t 'NLO',\
     fn u 1:($4/$5) w @NNLOline t 'NNLO',\
     fn u 1:($5/$5) w @N3LOline t 'N3LO'

set output
