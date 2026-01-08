set term postscript enhanced color 
set datafile fortran
set output 'bin_distribution.ps'

set ylabel 'N entries'
set xlabel 'picobarn'
file = 'index4_dist.top'

plot   file i 0 u (($1+$2)/2.):3:4 w e lw 2 ti 'index 4',\
       file i 1 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '68%',\
       file i 2 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '87%',\
       file i 3 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '92%',\
       file i 4 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '95%',\
       file i 5 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '97%',\
       file i 1 u (($1+$2)/2.):($3/2.):1:2 w xerrorbars lc rgb 'blue' lw 2 ti '1 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.2):((2.*$1)+(1.-2.)*(($1+$2)/2.)):((2.*$2)+(1.-2.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '2 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.4):((3.*$1)+(1.-3.)*(($1+$2)/2.)):((3.*$2)+(1.-3.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '3 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.7):((4.*$1)+(1.-4.)*(($1+$2)/2.)):((4.*$2)+(1.-4.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '4 sigma',\
       file i 1 u (($1+$2)/2.):($3/3.):((5.*$1)+(1.-5.)*(($1+$2)/2.)):((5.*$2)+(1.-5.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '5 sigma'


reset

set ylabel 'N entries'
set xlabel 'picobarn'
file = 'index7_dist.top'

plot   file i 0 u (($1+$2)/2.):3:4 w e lw 2 ti 'index 7',\
       file i 1 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '68%',\
       file i 2 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '87%',\
       file i 3 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '92%',\
       file i 4 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '95%',\
       file i 5 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '97%',\
       file i 1 u (($1+$2)/2.):($3/2.):1:2 w xerrorbars lc rgb 'blue' lw 2 ti '1 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.2):((2.*$1)+(1.-2.)*(($1+$2)/2.)):((2.*$2)+(1.-2.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '2 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.4):((3.*$1)+(1.-3.)*(($1+$2)/2.)):((3.*$2)+(1.-3.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '3 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.7):((4.*$1)+(1.-4.)*(($1+$2)/2.)):((4.*$2)+(1.-4.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '4 sigma',\
       file i 1 u (($1+$2)/2.):($3/3.):((5.*$1)+(1.-5.)*(($1+$2)/2.)):((5.*$2)+(1.-5.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '5 sigma'


reset

set ylabel 'N entries'
set xlabel 'picobarn'
file = 'index8_dist.top'

plot   file i 0 u (($1+$2)/2.):3:4 w e lw 2 ti 'index 8',\
       file i 1 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '68%',\
       file i 2 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '87%',\
       file i 3 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '92%',\
       file i 4 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '95%',\
       file i 5 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '97%',\
       file i 1 u (($1+$2)/2.):($3/2.):1:2 w xerrorbars lc rgb 'blue' lw 2 ti '1 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.2):((2.*$1)+(1.-2.)*(($1+$2)/2.)):((2.*$2)+(1.-2.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '2 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.4):((3.*$1)+(1.-3.)*(($1+$2)/2.)):((3.*$2)+(1.-3.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '3 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.7):((4.*$1)+(1.-4.)*(($1+$2)/2.)):((4.*$2)+(1.-4.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '4 sigma',\
       file i 1 u (($1+$2)/2.):($3/3.):((5.*$1)+(1.-5.)*(($1+$2)/2.)):((5.*$2)+(1.-5.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '5 sigma'


reset

set ylabel 'N entries'
set xlabel 'picobarn'
file = 'index9_dist.top'

plot   file i 0 u (($1+$2)/2.):3:4 w e lw 2 ti 'index 9',\
       file i 1 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '68%',\
       file i 2 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '87%',\
       file i 3 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '92%',\
       file i 4 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '95%',\
       file i 5 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '97%',\
       file i 1 u (($1+$2)/2.):($3/2.):1:2 w xerrorbars lc rgb 'blue' lw 2 ti '1 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.2):((2.*$1)+(1.-2.)*(($1+$2)/2.)):((2.*$2)+(1.-2.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '2 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.4):((3.*$1)+(1.-3.)*(($1+$2)/2.)):((3.*$2)+(1.-3.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '3 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.7):((4.*$1)+(1.-4.)*(($1+$2)/2.)):((4.*$2)+(1.-4.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '4 sigma',\
       file i 1 u (($1+$2)/2.):($3/3.):((5.*$1)+(1.-5.)*(($1+$2)/2.)):((5.*$2)+(1.-5.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '5 sigma'


reset

set ylabel 'N entries'
set xlabel 'picobarn'
file = 'index17_dist.top'

plot   file i 0 u (($1+$2)/2.):3:4 w e lw 2 ti 'index 17',\
       file i 1 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '68%',\
       file i 2 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '87%',\
       file i 3 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '92%',\
       file i 4 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '95%',\
       file i 5 u (($1+$2)/2.):3:1:2 w xerrorbars lc rgb 'green' lw 2 ti '97%',\
       file i 1 u (($1+$2)/2.):($3/2.):1:2 w xerrorbars lc rgb 'blue' lw 2 ti '1 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.2):((2.*$1)+(1.-2.)*(($1+$2)/2.)):((2.*$2)+(1.-2.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '2 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.4):((3.*$1)+(1.-3.)*(($1+$2)/2.)):((3.*$2)+(1.-3.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '3 sigma',\
       file i 1 u (($1+$2)/2.):($3/2.7):((4.*$1)+(1.-4.)*(($1+$2)/2.)):((4.*$2)+(1.-4.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '4 sigma',\
       file i 1 u (($1+$2)/2.):($3/3.):((5.*$1)+(1.-5.)*(($1+$2)/2.)):((5.*$2)+(1.-5.)*(($1+$2)/2.)) w xerrorbars lc rgb 'blue' lw 2 ti '5 sigma'


reset


set output

