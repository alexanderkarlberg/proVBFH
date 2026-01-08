set term postscript enhanced color
# set term wxt persist

set macros
calcerr = '(sqrt(($6/$5)**2+($14/$13)**2)*($5/$13))'


do for [type in '0-deriv 1-deriv 2-deriv 3-deriv 4-deriv 5-deriv 6-deriv'] {
set output 'checks-'.type.'.ps'
set grid
set logscale y
# set datafile fortran
set xlabel 'Y_{cut}'
set ylabel '{/Symbol s}_{wrong}/{/Symbol s}_{total}'
set xrange [0:5]
set yrange [1e-5:1]

dy=0.2

set title '11, 22 (NLO)'
plot '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 2 jets-11\") <(mergeidx.pl -f '.type.'.top \"sig tot 2 jets-11\")"' u ($3):($5/$13):(@calcerr) w errorl ti '2 jets, 11',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 2 jets-22\") <(mergeidx.pl -f '.type.'.top \"sig tot 2 jets-22\")"' u ($3):($5/$13):(@calcerr) w errorl ti '2 jets, 22',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-11\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-11\")"' u ($3):($5/$13):(@calcerr) w errorl ti '3 jets, 11',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-22\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-22\")"' u ($3):($5/$13):(@calcerr) w errorl ti '3 jets, 22',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-11\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-11\")"' u ($3):($5/$13):(@calcerr) w errorl ti '4 jets, 11',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-22\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-22\")"' u ($3):($5/$13):(@calcerr) w errorl lc rgb 'black' ti '4 jets, 22'

set title '12, 21 (NLO)'
plot '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 2 jets-12\") <(mergeidx.pl -f '.type.'.top \"sig tot 2 jets-12\")"' u ($3):($5/$13):(@calcerr) w errorl ti '2 jets, 12',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 2 jets-21\") <(mergeidx.pl -f '.type.'.top \"sig tot 2 jets-21\")"' u ($3):($5/$13):(@calcerr) w errorl ti '2 jets, 21',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-12\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-12\")"' u ($3):($5/$13):(@calcerr) w errorl ti '3 jets, 12',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-21\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-21\")"' u ($3):($5/$13):(@calcerr) w errorl ti '3 jets, 21',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-12\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-12\")"' u ($3):($5/$13):(@calcerr) w errorl ti '4 jets, 12',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-21\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-21\")"' u ($3):($5/$13):(@calcerr) w errorl lc rgb 'black' ti '4 jets, 21'

set yrange [1e-8:1.0]
set ylabel '{/Symbol s}_{wrong}'
set title '11, 22 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig wrong 2 jets-11"' u ($3):($5):($6) w errorl ti '2 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 2 jets-22"' u ($3):($5):($6) w errorl ti '2 jets, 22',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-11"' u ($3):($5):($6) w errorl ti '3 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-22"' u ($3):($5):($6) w errorl ti '3 jets, 22',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-11"' u ($3):($5):($6) w errorl ti '4 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-22"' u ($3):($5):($6) w errorl lc rgb 'black' ti '4 jets, 22'

set title '12, 21 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig wrong 2 jets-12"' u ($3):($5):($6) w errorl ti '2 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 2 jets-21"' u ($3):($5):($6) w errorl ti '2 jets, 21',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-12"' u ($3):($5):($6) w errorl ti '3 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-21"' u ($3):($5):($6) w errorl ti '3 jets, 21',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-12"' u ($3):($5):($6) w errorl ti '4 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-21"' u ($3):($5):($6) w errorl lc rgb 'black' ti '4 jets, 21'

set ylabel '{/Symbol s}_{total}'
set title '11, 22 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig tot 2 jets-11"' u ($3):($5):($6) w errorl ti '2 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig tot 2 jets-22"' u ($3):($5):($6) w errorl ti '2 jets, 22',\
     '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-11"' u ($3):($5):($6) w errorl ti '3 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-22"' u ($3):($5):($6) w errorl ti '3 jets, 22',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-11"' u ($3):($5):($6) w errorl ti '4 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-22"' u ($3):($5):($6) w errorl lc rgb 'black' ti '4 jets, 22'

set title '12, 21 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig tot 2 jets-12"' u ($3):($5):($6) w errorl ti '2 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig tot 2 jets-21"' u ($3):($5):($6) w errorl ti '2 jets, 21',\
     '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-12"' u ($3):($5):($6) w errorl ti '3 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-21"' u ($3):($5):($6) w errorl ti '3 jets, 21',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-12"' u ($3):($5):($6) w errorl ti '4 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-21"' u ($3):($5):($6) w errorl lc rgb 'black' ti '4 jets, 21'


set output
}
