set term postscript enhanced color
# set term wxt persist
set datafile fortran
set macros
calcerr = '(sqrt(($4/$3)**2+($8/$7)**2)*($3/$7))'


do for [type in 'nobug nobug_gen virt_bug virt_bug_gen_1d0 virt_bug_gen_005d0 sigcollremn_bug sigcollremn_bug_1d0 sigcollsoft_bug'] {
set output 'checks-'.type.'.ps'
set grid
set logscale y
# set datafile fortran
set xlabel 'Y_{cut}'
set ylabel '{/Symbol s}_{wrong}/{/Symbol s}_{total}'
set xrange [0:6]
set yrange [1e-4:1]

dy=0.2

set title '11, 22 (NLO)'
plot '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-11\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-11\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl ti '3 jets, 11',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-22\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-22\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl ti '3 jets, 22',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-11\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-11\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl ti '4 jets, 11',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-22\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-22\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl lc rgb 'black' ti '4 jets, 22'

set title '12, 21 (NLO)'
plot '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-12\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-12\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl ti '3 jets, 12',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 3 jets-21\") <(mergeidx.pl -f '.type.'.top \"sig tot 3 jets-21\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl ti '3 jets, 21',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-12\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-12\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl ti '4 jets, 12',\
     '< exec bash -c "paste <(mergeidx.pl -f '.type.'.top \"sig wrong 4 jets-21\") <(mergeidx.pl -f '.type.'.top \"sig tot 4 jets-21\")"' u (($1+$2)/2):($3/$7):(@calcerr) w errorl lc rgb 'black' ti '4 jets, 21'

set yrange [1e-8:1.0]
set ylabel '{/Symbol s}_{wrong}'
set title '11, 22 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-11"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-22"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 22',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-11"' u (($1+$2)/2):($3):($4) w errorl ti '4 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-22"' u (($1+$2)/2):($3):($4) w errorl lc rgb 'black' ti '4 jets, 22'

set title '12, 21 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-12"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 3 jets-21"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 21',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-12"' u (($1+$2)/2):($3):($4) w errorl ti '4 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig wrong 4 jets-21"' u (($1+$2)/2):($3):($4) w errorl lc rgb 'black' ti '4 jets, 21'

set ylabel '{/Symbol s}_{total}'
set title '11, 22 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-11"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-22"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 22',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-11"' u (($1+$2)/2):($3):($4) w errorl ti '4 jets, 11',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-22"' u (($1+$2)/2):($3):($4) w errorl lc rgb 'black' ti '4 jets, 22'

set title '12, 21 (NLO)'
plot '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-12"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig tot 3 jets-21"' u (($1+$2)/2):($3):($4) w errorl ti '3 jets, 21',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-12"' u (($1+$2)/2):($3):($4) w errorl ti '4 jets, 12',\
     '<mergeidx.pl -f '.type.'.top "sig tot 4 jets-21"' u (($1+$2)/2):($3):($4) w errorl lc rgb 'black' ti '4 jets, 21'


set output
}
