# gnuplot file

set term postscript enhanced colour dl 2 lw 2
filename='ptH.ps'
set output filename

set xlabel 'p_{tH} [GeV]'
set ylabel 'd{/Symbol s}/dp_{tH} [pb/GeV]'
set xrange [:400]
plot 'ptH-11.res' u 1:5 w l lt -1     t '11',\
     'ptH-22.res' u 1:5 w l lc 1 dt 2 t '22',\
     'ptH-HH.res' u 1:5 w l lc 1 dt 4 t 'HH',\
     'ptH-12.res' u 1:5 w l lc 2 dt 2 t '12',\
     'ptH-21.res' u 1:5 w l lc 2 dt 4 t '21',\
     'ptH-1H.res' u 1:5 w l lc 3 dt 2 t '1H',\
     'ptH-H1.res' u 1:5 w l lc 3 dt 4 t 'H1',\

set ylabel 'ratio to central scale'
plot '<paste ptH-11.res ptH-11.res' u 1:($10/$5) w l lt -1     t '11',\
     '<paste ptH-11.res ptH-22.res' u 1:($10/$5) w l lc 1 dt 2 t '22',\
     '<paste ptH-11.res ptH-HH.res' u 1:($10/$5) w l lc 1 dt 4 t 'HH',\
     '<paste ptH-11.res ptH-12.res' u 1:($10/$5) w l lc 2 dt 2 t '12',\
     '<paste ptH-11.res ptH-21.res' u 1:($10/$5) w l lc 2 dt 4 t '21',\
     '<paste ptH-11.res ptH-1H.res' u 1:($10/$5) w l lc 3 dt 2 t '1H',\
     '<paste ptH-11.res ptH-H1.res' u 1:($10/$5) w l lc 3 dt 4 t 'H1',\


# now investigate issue that total cross section doesn't agree with
# the integral over ptH
set auto x
dir='../../vbf_h3j_nnlo/paper_results/11_ptH/'
subfile=dir.'hjjj_subtraction_nlo.top'
incfile=dir.'hjj_inclusive_nnlo.top'
set ylabel '{/Symbol s}(below p_{tH})'

set title 'comparing exclusive part (subtraction): ptH integral to total X-sct'     
asympt=system('mergeidx.pl -f '.subfile." VBF | grep -v -e '#' -e '^$' | sed s/D/E/g | awk '{print $3}'")*1.0
plot '<mergeidx.pl -f '.subfile.' ptH | sed s/D/E/g | rebin.pl -2 -a -r -c 1' u 3:4 t 'ptH integral',\
     '<mergeidx.pl -f '.subfile.' ptj1 | sed s/D/E/g | rebin.pl -2 -a -r -c 1' u 3:4 t 'ptj1 integral',\
     '<mergeidx.pl -f '.subfile.' ptj2 | sed s/D/E/g | rebin.pl -2 -a -r -c 1' u 3:4 t 'ptj2 integral',\
     asympt t 'VBF total'

set title 'comparing inclusive part: ptH integral to total X-sct'     
asympt=system('mergeidx.pl -f '.incfile." VBF | grep -v -e '#' -e '^$' | sed s/D/E/g | awk '{print $3}'")*1.0
plot '<mergeidx.pl -f '.incfile.' ptH | sed s/D/E/g | rebin.pl -2 -a -r -c 1' u 3:4 t 'ptH integral',\
     '<mergeidx.pl -f '.incfile.' ptj1 | sed s/D/E/g | rebin.pl -2 -a -r -c 1' u 3:4 t 'ptj1 integral',\
     '<mergeidx.pl -f '.incfile.' ptj2 | sed s/D/E/g | rebin.pl -2 -a -r -c 1' u 3:4 t 'ptj2 integral',\
     asympt t 'VBF total'


set output
system("pstopdf ".filename)     
