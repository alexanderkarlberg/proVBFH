set term pdfcairo enhanced color transparent font "Times, 26" size 20cm,15cm
# set term wxt persist

fn = 'n3lo_hist.dat'


colLO   = '#FFCC11'
colNLO  = '#0099CC'
colNNLO = '#37BC61'
colN3LO = '#FF4500'

# colLO   = '#FFF8DC'
# colNLO  = '#BFEFFF'
# colNNLO = '#37BC61'
# colN3LO = '#FF4500'

set output 'n3lo-total.pdf'

# plot ptH
reset
set grid
set xtics format ""
set ytics 0.001,10,100
set ylabel '{/Symbol s} [pb]'
set label 2000 'NNPDF30\_nnlo\_as\_0118' at graph 0.02,0.07
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.02,0.15
set label 2002 'LHC 13 TeV' at graph 0.02,0.22
set label 1000 '{/*4 PRELIMINARY' at graph 0.45,0.50 center rotate by 25 tc rgb '#d0d0d0' back

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
set grid front
set label 1000 '{/*4 PRELIMINARY' at graph 0.45,0.50 center rotate by 25 tc rgb '#d0d0d0' 
set label 2000 'NNPDF30\_nnlo\_as\_0118' at graph 0.02,0.07
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.02,0.15
set label 2002 'LHC 13 TeV' at graph 0.02,0.22
set key maxrows 2
set key top right
set key width 3
set format y "10^{%L}"
set ytics 0.001,10,100
set xtics format ""
set ylabel 'd{/Symbol s}/dp_{t,H} [pb]'

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
set yrange [1e-4:1e-1]
set label 1000 back
plot '<mergeidx.pl -f '.fn.' LO.ptH'   u (($1+$2)/2):3 w lines  lc rgb colLO lw 2.0 not,\
     '' u (($1+$2)/2):4:5  w filledcurves lc rgb colLO lw 2.0 fs transparent solid 0.5 ti 'LO',\
     '<mergeidx.pl -f '.fn.' NLO.ptH'  u (($1+$2)/2):3 w lines lc rgb colNLO lw 2.0 not,\
     '' u (($1+$2)/2.):4:5  w filledcurves lc rgb colNLO lw 2.0 fs transparent solid 0.5 ti 'NLO',\
     '<mergeidx.pl -f '.fn.' NNLO.ptH' u (($1+$2)/2):3 w lines lc rgb colNNLO lw 2.0 not,\
     '' u (($1+$2)/2):4  w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2):5 w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2.):4:5 w filledcurves lc rgb colNNLO lw 2.0 fs transparent solid 0.5 ti 'NNLO',\
     '<mergeidx.pl -f '.fn.' N3LO.ptH' u (($1+$2)/2):3 w lines lc rgb colN3LO lw 2.0 not,\
     '' u (($1+$2)/2):4  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):5 w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2.):4:5 w filledcurves lc rgb colN3LO lw 2.0 fs transparent solid 0.5 ti 'N3LO'


set tmargin at screen 0.47
set bmargin at screen 0.115
set nologscale y
set title ''
set ytics 0.9,0.02,1.2
set format x
set format y
unset label
set ylabel 'ratio to N3LO'
set xlabel 'p_{t,H} [GeV]'
set yrange [0.98:1.05]
set xrange [0:300]
plot '<mergeidx.pl -f '.fn.' LO.ptH N3LO.ptH'    u (($1+$2)/2):($3/$8) w lines  lc rgb colLO lw 2.0 not,\
     '' u (($1+$2)/2):($4/$8):($5/$8)  w filledcurves lc rgb colLO lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl -f '.fn.' NLO.ptH N3LO.ptH'  u (($1+$2)/2):($3/$8) w lines lc rgb colNLO lw 2.0 not,\
     '' u (($1+$2)/2.):($4/$8):($5/$8)  w filledcurves lc rgb colNLO lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl -f '.fn.' NNLO.ptH N3LO.ptH' u (($1+$2)/2):($3/$8) w lines lc rgb colNNLO lw 2.0 not,\
     '' u (($1+$2)/2):($4/$8)  w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2):($5/$8) w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2.):($4/$8):($5/$8) w filledcurves lc rgb colNNLO lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl -f '.fn.' N3LO.ptH'  u (($1+$2)/2):3 w lines lc rgb colN3LO lw 2.0 not,\
     '' u (($1+$2)/2):($4/$3)  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):($5/$3) w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2.):($4/$3):($5/$3) w filledcurves lc rgb colN3LO lw 2.0 fs transparent solid 0.5 not

unset multiplot

# YH plot
set grid front
set label 1000 '{/*4 PRELIMINARY' at graph 0.45,0.50 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'NNPDF30\_nnlo\_as\_0118' at graph 0.02,0.07
set label 2001 'Q/2 < {/Symbol m}_R , {/Symbol m}_F < 2 Q' at graph 0.02,0.15
set label 2002 'LHC 13 TeV' at graph 0.02,0.22
set key maxrows 2
set key top right
set key width 3
set format y "10^{%L}"
set xtics format ""
set xlabel ''
set ylabel 'd{/Symbol s}/dy_H [pb]'
set multiplot
set origin 0.0,0.3
set size 1.0,1.0
iindex = 7
set lmargin at screen 0.20
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.47
set logscale y
set ytics 0.001,10,100
set yrange [1e-4:1e1]
set xrange [0:4.5]
plot '<mergeidx.pl -f '.fn.' LO.yH'   u (($1+$2)/2):3 w lines  lc rgb colLO lw 2.0 not,\
     '' u (($1+$2)/2):4:5  w filledcurves lc rgb colLO lw 2.0 fs transparent solid 0.5 ti 'LO',\
     '<mergeidx.pl -f '.fn.' NLO.yH'  u (($1+$2)/2):3 w lines lc rgb colNLO lw 2.0 not,\
     '' u (($1+$2)/2.):4:5  w filledcurves lc rgb colNLO lw 2.0 fs transparent solid 0.5 ti 'NLO',\
     '<mergeidx.pl -f '.fn.' NNLO.yH' u (($1+$2)/2):3 w lines lc rgb colNNLO lw 2.0 not,\
     '' u (($1+$2)/2):4  w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2):5 w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2.):4:5 w filledcurves lc rgb colNNLO lw 2.0 fs transparent solid 0.5 ti 'NNLO',\
     '<mergeidx.pl -f '.fn.' N3LO.yH' u (($1+$2)/2):3 w lines lc rgb colN3LO lw 2.0 not,\
     '' u (($1+$2)/2):4  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):5 w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2.):4:5 w filledcurves lc rgb colN3LO lw 2.0 fs transparent solid 0.5 ti 'N3LO'

set tmargin at screen 0.47
set bmargin at screen 0.115
set nologscale y
set title ''
set ytics 0.9,0.02,1.2
set format x
set format y
set ylabel 'ratio to N3LO'
unset label
set xlabel '|y_{H}| [GeV]'
set yrange [0.98:1.05]
set xrange [0:4.5]
plot '<mergeidx.pl -f '.fn.' LO.yH N3LO.yH'    u (($1+$2)/2):($3/$8) w lines  lc rgb colLO lw 2.0 not,\
     '' u (($1+$2)/2):($4/$8):($5/$8)  w filledcurves lc rgb colLO lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl -f '.fn.' NLO.yH N3LO.yH'  u (($1+$2)/2):($3/$8) w lines lc rgb colNLO lw 2.0 not,\
     '' u (($1+$2)/2.):($4/$8):($5/$8)  w filledcurves lc rgb colNLO lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl -f '.fn.' NNLO.yH N3LO.yH' u (($1+$2)/2):($3/$8) w lines lc rgb colNNLO lw 2.0 not,\
     '' u (($1+$2)/2):($4/$8)  w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2):($5/$8) w lines lc rgb colNNLO lw 0.5 not,\
     '' u (($1+$2)/2.):($4/$8):($5/$8) w filledcurves lc rgb colNNLO lw 2.0 fs transparent solid 0.5 not,\
     '<mergeidx.pl -f '.fn.' N3LO.yH'  u (($1+$2)/2):($3/$3) w lines lc rgb colN3LO lw 2.0 not,\
     '' u (($1+$2)/2):($4/$3)  w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2):($5/$3) w lines lc rgb colN3LO lw 0.5 not,\
     '' u (($1+$2)/2.):($4/$3):($5/$3) w filledcurves lc rgb colN3LO lw 2.0 fs transparent solid 0.5 not

unset multiplot

set output
