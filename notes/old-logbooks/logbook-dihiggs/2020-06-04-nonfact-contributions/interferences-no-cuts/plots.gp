set terminal pdf enhanced font "Latin Modern Roman,20" size 18cm,25cm
set style fill border
set key noautotitles
#set tmargin 2
#set bmargin 0.5
rsize = 0.05
tsize = 0.05

set output "contributions-no-cuts.pdf"
set datafile fortran

full='full-sum.top'
box='box-sum.top'
tri='tri-only.top'
box_t='box-t-only.top'
box_u='box-u-only.top'
full_interference='interference-full.top'
box_t_u_interference='box-t-u-interference.top'
tri_box_interference='tri-box-interference.top'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=4

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=5

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=6

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=7

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=8

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=9

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=10

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=11

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=12

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=13

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=14

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=15

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=16

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=17

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=18

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=19

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=20

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=21

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=22

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'

reset
set size 1,0.97
set origin 0,0.0
set key vertical right top invert
set key spacing 1.1 width 2 height 1

set style fill transparent solid 0.4 noborder

set mxtics
set mytics
set grid

set xrange [*:*]
set yrange [*:*]
set grid
set style fill transparent pattern 6 border
iindex=23

pl full i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Full NNLO nonfact',\
   box i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box only',\
   box_t i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t only',\
   box_u i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box u only',\
   tri i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle only',\
   full_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'All interference',\
   tri_box_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Triangle/box interference',\
   box_t_u_interference i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'Box t/u interference'


set output

