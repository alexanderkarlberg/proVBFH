set term postscript enhanced color
set datafile fortran
# set term wxt persist
set output 'scales.ps'

# output from altern_method = 2
fn_max = 'maxQ.top'
fn_min = 'minQ.top'
fn_sqrt = 'sqrtQ1Q2.top'

# plot ptH
reset
set grid
set xtics format ""
set xrange [0:400]
set format x
#set yrange [1e-6:5e-2]
set title 'Q1,Q2 comparisons'
set xlabel "p_{T,H} [GeV]"
set ylabel "[Gev]"
iindex = 0

mH=125.0
mW=80.385
mZ=91.1876
mWH=sqrt(mH*mW)
Qmaxstr='((mH/2)**4 + (ptH)**4)**(1.0/4)'
Qminstr= 'mWH/2 + (mZ - mWH/2)*(ptH**2)/(ptH**2+mH**2)'

# new choices
Qmaxstr='sqrt((mH/2)**2 + ptH**2)'
Qminstr='mH/2'
# an alternative "single" choice
Qmidstr='sqrt((mH/2)**2 + mH*ptH/2)'

eval 'Qmax(ptH)='.Qmaxstr
eval 'Qmin(ptH)='.Qminstr
eval 'Qmid(ptH)='.Qmidstr

set label 1 at graph 0.05,0.95 'Qmax(ptH)='.Qmaxstr
set label 2 at graph 0.05,0.88 'Qmin(ptH)='.Qminstr
set label 3 at graph 0.05,0.81 'Qmid(ptH)='.Qmidstr

plot fn_max i iindex u (($1+$2)/2):3 linewidth 1.0 ti 'max{Q1,Q2}',\
     fn_min i iindex u (($1+$2)/2):3 linewidth 1.0 ti 'min{Q1,Q2}',\
     fn_sqrt i iindex u (($1+$2)/2):3 linewidth 1.0 ti 'sqrt(Q1*Q2)',\
     Qmid(x) lw 3 lc rgb 'red',\
     Qmin(x) lw 4,\
     Qmax(x) lw 4,\
     sqrt(Qmin(x)*Qmax(x)) lw 4 lc 0

#     exp(-20.1/x)*(x*80.385)**0.5 + exp(-x/31.25)*(125/2+x) lw 3 lc rgb 'red',\

# plot ptH
#reset
set grid
set xtics format ""
set xrange [400:2000]
set log x
#set log y
set format x
#set yrange [1e-6:5e-2]
set title 'Q1,Q2 comparisons'
set xlabel "p_{T,H} [GeV]"
set ylabel "[Gev]"
iindex = 0

mH=125.0
mW=81.285
mZ=90.1876
mWH=sqrt(mH*mW)

plot fn_max i iindex u (($1+$2)/2):3 linewidth 1.0 ti 'max{Q1,Q2}',\
     fn_min i iindex u (($1+$2)/2):3 linewidth 1.0 ti 'min{Q1,Q2}',\
     fn_sqrt i iindex u (($1+$2)/2):3 linewidth 1.0 ti 'sqrt(Q1*Q2)',\
     Qmid(x) lw 3 lc rgb 'red',\
     Qmin(x) lw 2,\
     Qmax(x) lw 2,\
     sqrt(Qmin(x)*Qmax(x)) lw 2 lc 0

#     exp(-20.1/x)*(x*80.385)**0.5 + exp(-x/31.25)*(125/2+x) lw 5 lc rgb 'red',\

# get the same plot on a log scale
set log
set xrange [10:]
set yrange [40:3000]
replot

set output


# gnuplot file

#set xrange [0:400]

#mH=125.0
#mW=81.285
#mZ=90.1876
#mWH=sqrt(mH*mW)
#Qmax(ptH)=((mH/2)**4 + (ptH)**4)**(1.0/4)
#Qmin(ptH)= mWH/2 + (mZ - mWH/2)*(ptH**2)/(ptH**2+mH**2)

#set yrange [0:400]
#set grid

#plot Qmin(x), Qmax(x), sqrt(Qmin(x)*Qmax(x))
