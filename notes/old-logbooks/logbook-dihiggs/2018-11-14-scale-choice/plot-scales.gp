set term pdfcairo enhanced color transparent font "Times, 26" size 22cm,17cm

set output 'test-scales.pdf'

set xlabel "p_{t,HH} [GeV]"
set ylabel "[Gev]"
iindex = 0

mH=125.0
mW=80.385
mZ=91.1876

# new choices
Qmaxstr='sqrt((mH/2)**2 + ptHH**2)'
Qminstr='mH/2'
# an alternative "single" choice
Qmidstr='sqrt((mH/2)**2 + mH*ptHH/2)'
Qmidbisstr='sqrt((mH)**2 + mH*ptHH)'

eval 'Qmax(ptHH)='.Qmaxstr
eval 'Qmin(ptHH)='.Qminstr
eval 'Qmid(ptHH)='.Qmidstr
eval 'Qmidbis(ptHH)='.Qmidbisstr

set label 2 at graph 0.98,0.10 '{/Symbol m}@^2_{max}(p_{t,HH}) = m@^2_h/4 + p@^2_{t,HH}' right
set label 1 at graph 0.98,0.04 '{/Symbol m}@^2_0(p_{t,HH}) = m@^2_h/4 + m_h p_{t,HH}/2' right
mrg='<mergeidx.pl -f test-scales.top '
set key top left
set yrange [0:350]
plot mrg.'sqrt.Q1.Q2 norm' u (0.5*($1+$2)):($3/$7) w l lc 1 lw 2 t 'sqrt(Q_1 Q_2)',\
     mrg.'max.Q1.Q2 norm'  u (0.5*($1+$2)):($3/$7) w l lc 3 lw 2 dt 2 t 'max(Q_1,Q_2)',\
     mrg.'min.Q1.Q2 norm'  u (0.5*($1+$2)):($3/$7) w l lc 3 lw 2 dt (4,6) t 'min(Q_1,Q_2)',\
     mrg.'ptHH norm'       u (0.5*($1+$2)):($3/$7) w l lc 4 lw 2 dt 2 t 'p_{t,HH}',\
     Qmin(x) lc 4 lw 2 dt (4,6) t 'm_h/2',\
     Qmid(x) w l lc 7 lw 2 t '{/Symbol m}_0(p_{t,HH})'

set yrange [0:400]
set key maxrow 5
set label 2 at graph 0.98,0.10 '{/Symbol m}@^2_{max} = m@^2_h/4 + p@^2_{t,HH}' right
set label 1 at graph 0.98,0.04 '{/Symbol m}@^2_0 = m@^2_h/4 + m_h p_{t,HH}/2' right
set label 4 at graph 0.63,0.10 '{/Symbol m}@^2_{alt} = m@^2_h + m_h p_{t,HH}' right
set label 3 at graph 0.63,0.04 '{/Symbol m}@^2_{dyn} = m@^2_{HH}/4 + m_{HH} p_{t,HH}/2' right
replot Qmax(x) lc 2 lw 2 dt (4,6) t '{/Symbol m}_{max}',\
       Qmidbis(x) lw 2 dt 2 t '{/Symbol m}_{alt}',\
       mrg.'muDyn norm' u (0.5*($1+$2)):($3/$7) w l lc 7 lw 2 dt 2 t '{/Symbol m}_{dyn}(m_{HH},p_{t,HH})'
