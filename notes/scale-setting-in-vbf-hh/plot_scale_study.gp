set term pdfcairo enhanced color transparent font "Times, 34" size 30cm,22cm

fn = 'pwg-LO-0001.top'

set output 'scale_study.pdf'
set macros
mH = 125
Qmin(x) = mH/2 
Qmax(x) = sqrt((mH/2)**2 + x**2)
mu(x) = sqrt(Qmin(x) * Qmax(x))
LOline   = "lines lc rgb '#f9B919' lw 3.0"
LOlinethin   = "lines lc rgb '#fccc55' lw 0.5"
NLOline  = "lines lc rgb '#38B0DE' lw 3.0"
NLOlinethin  = "lines lc rgb '#05AEEFFF' lw 0.5"
NNLOline = "lines lc rgb '#37BC61' lw 3.0"
NNLOlinethin = "lines lc rgb '#37BC61' lw 0.5"
N3LOline = "lines lc rgb '#FF4500' lw 3.0"

LOfill   = "filledcurves lc rgb '#fccc55' lw 2.0 fs transparent pattern 4"
#LOfill   = "filledcurves lc rgb '#FFF0C4' lw 2.0 fs transparent solid 0.8"
NLOfill  = "filledcurves lc rgb '#05AEEFFF' lw 2.0 fs transparent pattern 5"
# NLOfill  = "filledcurves lc rgb '#05AEEFFF' lw 2.0 fs transparent solid 0.85"
NNLOfill = "filledcurves lc rgb '#37BC61' lw 2.0 fs transparent pattern 1"
N3LOfill = "filledcurves lc rgb '#FF4500' lw 2.0 fs transparent solid 0.5"

# plot cross section
reset
#set grid front #mytics mxtics ytics xtics
#set logscale y
#set xtics format ""
set ylabel '<μ> [GeV]'
set xlabel 'p_{t,HH} [GeV]'
# set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
#set label 2004 'SM' at graph 0.51,0.94
set label 2003 'proVBFHH\@LO v2.1.0' at graph 0.63,0.94
set label 2002 'LHC 14 TeV' at graph 0.63,0.88
set label 2001 'PDF4LHC21\_40' at graph 0.63,0.82
#set label 2000 'Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.03,0.06
#set key maxrows 2
set key top left
set key width 1
set key spacing 1.35
#set ytics 0.1,10,1000000
set xtics 0,50,400
set mxtics 5
set yrange [40:450]
set xrange [*:*]
set grid
#set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'black' lw 3 dt 2

plot '<mergeidx.pl ptHH sqrtQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 1 lc rgb 'red' title '{/Symbol=\326}Q_1Q_2',\
     '<mergeidx.pl ptHH minQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 2 lc rgb 'red' title 'min\{Q_1,Q_2\}',\
     '<mergeidx.pl ptHH maxQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 3 lc rgb 'red' title 'max\{Q_1,Q_2\}',\
     Qmin(x) lw 6 dt 2 lc rgb 'black' t 'm_H/2',\
     Qmax(x) lw 6 dt 3 lc rgb 'black' t '{/Symbol=\326}(m_H/2)^2 + p@_{t,HH}^2',\
     mu(x) lw 6 dt 1 lc rgb 'black' t '{/Symbol=\326}m_H/2{/Symbol=\326}(m_H/2)^2 + p@_{t,HH}^2',\
     '<mergeidx.pl ptHH sqrtptj1ptj2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 2 lc rgb '#f9B919' title '{/Symbol=\326}p_{t,j_1}p_{t,j_2}',\
     '<mergeidx.pl ptHH sqrtMTH1MTH2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3*0.5) w l lw 6 dt 2 lc rgb '#38B0DE' title '({/Symbol=\326}M_{t,H_1}M_{t,H_2})/2'

reset
#set grid front #mytics mxtics ytics xtics
#set logscale y
#set xtics format ""
set ylabel '<μ> [GeV]'
set xlabel 'p_{t,HH} [GeV]'
# set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
#set label 2004 'SM' at graph 0.51,0.94
set label 2003 'proVBFHH\@LO v2.1.0' at graph 0.63,0.94
set label 2002 'LHC 14 TeV' at graph 0.63,0.88
set label 2001 'PDF4LHC21\_40' at graph 0.63,0.82
#set label 2000 'Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.03,0.06
#set key maxrows 2
set key top left
set key width 1
set key spacing 1.35
#set ytics 0.1,10,1000000
set xtics 0,50,400
set mxtics 5
set yrange [40:180]
set xrange [*:150]
set grid
#set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'black' lw 3 dt 2

plot '<mergeidx.pl ptHH sqrtQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 1 lc rgb 'red' title '{/Symbol=\326}Q_1Q_2',\
     '<mergeidx.pl ptHH minQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 2 lc rgb 'red' title 'min\{Q_1,Q_2\}',\
     '<mergeidx.pl ptHH maxQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 3 lc rgb 'red' title 'max\{Q_1,Q_2\}',\
     Qmin(x) lw 6 dt 2 lc rgb 'black' t 'm_H/2',\
     Qmax(x) lw 6 dt 3 lc rgb 'black' t '{/Symbol=\326}(m_H/2)^2 + p@_{t,HH}^2',\
     mu(x) lw 6 dt 1 lc rgb 'black' t '{/Symbol=\326}m_H/2{/Symbol=\326}(m_H/2)^2 + p@_{t,HH}^2',\
     '<mergeidx.pl ptHH sqrtptj1ptj2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 2 lc rgb '#f9B919' title '{/Symbol=\326}p_{t,j_1}p_{t,j_2}',\
     '<mergeidx.pl ptHH sqrtMTH1MTH2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3*0.5) w l lw 6 dt 2 lc rgb '#38B0DE' title '({/Symbol=\326}M_{t,H_1}M_{t,H_2})/2'

reset
#set grid front #mytics mxtics ytics xtics
#set logscale y
#set xtics format ""
set ylabel '<μ> [GeV]'
set xlabel 'p_{t,HH} [GeV]'
# set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
#set label 2004 'SM' at graph 0.51,0.94
set label 2003 'proVBFHH\@LO v2.1.0' at graph 0.63,0.94
set label 2002 'LHC 14 TeV' at graph 0.63,0.88
set label 2001 'PDF4LHC21\_40' at graph 0.63,0.82
#set label 2000 'Q_i/2 < {/Symbol m}_{R,i} , {/Symbol m}_{F,i} < 2 Q_i' at graph 0.03,0.06
#set key maxrows 2
set key top left
set key width 1
set key spacing 1.35
#set ytics 0.1,10,1000000
set xtics 0,50,400
set mxtics 5
set yrange [40:250]
set xrange [*:*]
set grid
#set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'black' lw 3 dt 2

plot '<mergeidx.pl ptHH sqrtQ1Q2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 1 lc rgb 'red' title '{/Symbol=\326}Q_1Q_2',\
     mu(x) lw 6 dt 1 lc rgb 'black' t '{/Symbol=\326}m_H/2{/Symbol=\326}(m_H/2)^2 + p@_{t,HH}^2',\
     '<mergeidx.pl ptHH sqrtptj1ptj2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3) w l lw 6 dt 2 lc rgb '#f9B919' title '{/Symbol=\326}p_{t,j_1}p_{t,j_2}',\
     '<mergeidx.pl ptHH sqrtMTH1MTH2 -f pwg-LO-0001.top' u (($1+$2)/2.):($7/$3*0.5) w l lw 6 dt 2 lc rgb '#38B0DE' title '({/Symbol=\326}M_{t,H_1}M_{t,H_2})/2'

set output
