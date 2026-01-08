# gnuplot file


reset
set sty data linespoints

set term postscript enhanced colour solid lw 1

set macros

badfactor=0.995
LOWp='($12-$2)/$2'
NLOWp='($12+$14-$2-$4)/$2'
NNLOWp='($12+$14+$16-$2-$4-$6)/$2'
NNLObadWp=sprintf('($12+$14+$16-$2-$4-%f*$6)/$2',badfactor)
LOWm='($13-$3)/$3'
NLOWm='($13+$15-$3-$5)/$3'
NNLOWm='($13+$15+$17-$3-$5-$7)/$3'
NNLObadWm=sprintf('($13+$15+$17-$3-$5-%f*$7)/$3',badfactor)
LOZ='($18-$8)/$8'
NLOZ='($18+$19-$8-$9)/$8'
NNLOZ='($18+$19+$20-$8-$9-$10)/$8'
NNLObadZ=sprintf('($18+$19+$20-$8-$9-%f*$10)/$8',badfactor)

do for [structf in '1 2 3']{
filename='check-scale-dep-f'.structf.'Wp-originalcode-prec1e-6.ps'
set output filename

# loop over scale choices in format
# ab_cd => xmuR1=a, xmuF1=b, xmuR2=c, xmuF2=d
scalechoices='11_14 11_41 11_43 14_41 14_43 41_43'

do for [scale in scalechoices] {
xR1 = scale[1:1]
xF1 = scale[2:2]
xR2 = scale[4:4]
xF2 = scale[5:5]
print 'xmuR1 = '.xR1.', xmuF1 = '.xF1.', xmuR2 = '.xR2.', xmuF2 = '.xF2
base(r,f)='original-F'.structf.'WpmZ_mur'.r.'.0_muf'.f.'.0_scQ.dat'
filepair(as)='<paste tmp-as'.as.'-prec1e-6/'.base(xR1,xF1).' tmp-as'.as.'-prec1e-6/'.base(xR2,xF2)
set log x

set xrange [0.01:0.8]
set xlabel 'x'
set ylabel '(scale-dependence/LO-result) * 1/{/Symbol a}^{p+1} [arb. units]'


set label 1 at graph 0.05,0.05
set label 2 at graph 0.05,0.10 '{/*0.7 '.system('echo "'.filepair('0.1024').'" | sed "s/_/\\\_/g"')
set label 3 at graph 0.05,0.15 'As {/Symbol a} -> 0, curves should tend to overlap'
set label 4 at graph 0.05,0.20 '(xmuR, xmuF) = ('.xR1.','.xF1.'), ('.xR2.','.xF2.')'

pow=1
set title 'Check at LO'
set label 1 LOWp
plot filepair('0.1024') u 1:(  1**pow*@LOWp),\
     filepair('0.0512') u 1:(  2**pow*@LOWp),\
     filepair('0.0256') u 1:(  4**pow*@LOWp),\
     filepair('0.0128') u 1:(  8**pow*@LOWp),\
     filepair('0.0064') u 1:( 16**pow*@LOWp),\
     filepair('0.0032') u 1:( 32**pow*@LOWp),\
     filepair('0.0016') u 1:( 64**pow*@LOWp),\
     filepair('0.0008') u 1:(128**pow*@LOWp)


pow=2
set title 'Check at NLO'
set label 1 NLOWp
plot filepair('0.1024') u 1:(  1**pow*@NLOWp),\
     filepair('0.0512') u 1:(  2**pow*@NLOWp),\
     filepair('0.0256') u 1:(  4**pow*@NLOWp),\
     filepair('0.0128') u 1:(  8**pow*@NLOWp),\
     filepair('0.0064') u 1:( 16**pow*@NLOWp),\
     filepair('0.0032') u 1:( 32**pow*@NLOWp),\
     filepair('0.0016') u 1:( 64**pow*@NLOWp),\
     filepair('0.0008') u 1:(128**pow*@NLOWp)


pow=3
set title 'Check at NNLO'
set label 1 NNLOWp
plot filepair('0.1024') u 1:(  1**pow*@NNLOWp),\
     filepair('0.0512') u 1:(  2**pow*@NNLOWp),\
     filepair('0.0256') u 1:(  4**pow*@NNLOWp),\
     filepair('0.0128') u 1:(  8**pow*@NNLOWp),\
     filepair('0.0064') u 1:( 16**pow*@NNLOWp),\
     filepair('0.0032') u 1:( 32**pow*@NNLOWp),\
     filepair('0.0016') u 1:( 64**pow*@NNLOWp),\
     filepair('0.0008') u 1:(128**pow*@NNLOWp)

set title sprintf('Check of artificial bad NNLO (with NNLO piece scaled by %f)',badfactor)
set label 1 NNLObadWp
plot filepair('0.1024') u 1:(  1**pow*@NNLObadWp),\
     filepair('0.0512') u 1:(  2**pow*@NNLObadWp),\
     filepair('0.0256') u 1:(  4**pow*@NNLObadWp),\
     filepair('0.0128') u 1:(  8**pow*@NNLObadWp),\
     filepair('0.0064') u 1:( 16**pow*@NNLObadWp),\
     filepair('0.0032') u 1:( 32**pow*@NNLObadWp),\
     filepair('0.0016') u 1:( 64**pow*@NNLObadWp),\
     filepair('0.0008') u 1:(128**pow*@NNLObadWp)
}
set output
}

do for [structf in '1 2 3']{
filename='check-scale-dep-f'.structf.'Wm-originalcode-prec1e-6.ps'
set output filename

# loop over scale choices in format
# ab_cd => xmuR1=a, xmuF1=b, xmuR2=c, xmuF2=d
scalechoices='11_14 11_41 11_43 14_41 14_43 41_43'

do for [scale in scalechoices] {
xR1 = scale[1:1]
xF1 = scale[2:2]
xR2 = scale[4:4]
xF2 = scale[5:5]
print 'xmuR1 = '.xR1.', xmuF1 = '.xF1.', xmuR2 = '.xR2.', xmuF2 = '.xF2
base(r,f)='original-F'.structf.'WpmZ_mur'.r.'.0_muf'.f.'.0_scQ.dat'
filepair(as)='<paste tmp-as'.as.'-prec1e-6/'.base(xR1,xF1).' tmp-as'.as.'-prec1e-6/'.base(xR2,xF2)
set log x

set xrange [0.01:0.8]
set xlabel 'x'
set ylabel '(scale-dependence/LO-result) * 1/{/Symbol a}^{p+1} [arb. units]'


set label 1 at graph 0.05,0.05
set label 2 at graph 0.05,0.10 '{/*0.7 '.system('echo "'.filepair('0.1024').'" | sed "s/_/\\\_/g"')
set label 3 at graph 0.05,0.15 'As {/Symbol a} -> 0, curves should tend to overlap'
set label 4 at graph 0.05,0.20 '(xmuR, xmuF) = ('.xR1.','.xF1.'), ('.xR2.','.xF2.')'

pow=1
set title 'Check at LO'
set label 1 LOWm
plot filepair('0.1024') u 1:(  1**pow*@LOWm),\
     filepair('0.0512') u 1:(  2**pow*@LOWm),\
     filepair('0.0256') u 1:(  4**pow*@LOWm),\
     filepair('0.0128') u 1:(  8**pow*@LOWm),\
     filepair('0.0064') u 1:( 16**pow*@LOWm),\
     filepair('0.0032') u 1:( 32**pow*@LOWm),\
     filepair('0.0016') u 1:( 64**pow*@LOWm),\
     filepair('0.0008') u 1:(128**pow*@LOWm)


pow=2
set title 'Check at NLO'
set label 1 NLOWm
plot filepair('0.1024') u 1:(  1**pow*@NLOWm),\
     filepair('0.0512') u 1:(  2**pow*@NLOWm),\
     filepair('0.0256') u 1:(  4**pow*@NLOWm),\
     filepair('0.0128') u 1:(  8**pow*@NLOWm),\
     filepair('0.0064') u 1:( 16**pow*@NLOWm),\
     filepair('0.0032') u 1:( 32**pow*@NLOWm),\
     filepair('0.0016') u 1:( 64**pow*@NLOWm),\
     filepair('0.0008') u 1:(128**pow*@NLOWm)


pow=3
set title 'Check at NNLO'
set label 1 NNLOWm
plot filepair('0.1024') u 1:(  1**pow*@NNLOWm),\
     filepair('0.0512') u 1:(  2**pow*@NNLOWm),\
     filepair('0.0256') u 1:(  4**pow*@NNLOWm),\
     filepair('0.0128') u 1:(  8**pow*@NNLOWm),\
     filepair('0.0064') u 1:( 16**pow*@NNLOWm),\
     filepair('0.0032') u 1:( 32**pow*@NNLOWm),\
     filepair('0.0016') u 1:( 64**pow*@NNLOWm),\
     filepair('0.0008') u 1:(128**pow*@NNLOWm)

set title sprintf('Check of artificial bad NNLO (with NNLO piece scaled by %f)',badfactor)
set label 1 NNLObadWm
plot filepair('0.1024') u 1:(  1**pow*@NNLObadWm),\
     filepair('0.0512') u 1:(  2**pow*@NNLObadWm),\
     filepair('0.0256') u 1:(  4**pow*@NNLObadWm),\
     filepair('0.0128') u 1:(  8**pow*@NNLObadWm),\
     filepair('0.0064') u 1:( 16**pow*@NNLObadWm),\
     filepair('0.0032') u 1:( 32**pow*@NNLObadWm),\
     filepair('0.0016') u 1:( 64**pow*@NNLObadWm),\
     filepair('0.0008') u 1:(128**pow*@NNLObadWm)
}
set output
}

do for [structf in '1 2 3']{
filename='check-scale-dep-f'.structf.'Z-originalcode-prec1e-6.ps'
set output filename

# loop over scale choices in format
# ab_cd => xmuR1=a, xmuF1=b, xmuR2=c, xmuF2=d
scalechoices='11_14 11_41 11_43 14_41 14_43 41_43'

do for [scale in scalechoices] {
xR1 = scale[1:1]
xF1 = scale[2:2]
xR2 = scale[4:4]
xF2 = scale[5:5]
print 'xmuR1 = '.xR1.', xmuF1 = '.xF1.', xmuR2 = '.xR2.', xmuF2 = '.xF2
base(r,f)='original-F'.structf.'WpmZ_mur'.r.'.0_muf'.f.'.0_scQ.dat'
filepair(as)='<paste tmp-as'.as.'-prec1e-6/'.base(xR1,xF1).' tmp-as'.as.'-prec1e-6/'.base(xR2,xF2)
set log x

set xrange [0.01:0.8]
set xlabel 'x'
set ylabel '(scale-dependence/LO-result) * 1/{/Symbol a}^{p+1} [arb. units]'


set label 1 at graph 0.05,0.05
set label 2 at graph 0.05,0.10 '{/*0.7 '.system('echo "'.filepair('0.1024').'" | sed "s/_/\\\_/g"')
set label 3 at graph 0.05,0.15 'As {/Symbol a} -> 0, curves should tend to overlap'
set label 4 at graph 0.05,0.20 '(xmuR, xmuF) = ('.xR1.','.xF1.'), ('.xR2.','.xF2.')'

pow=1
set title 'Check at LO'
set label 1 LOZ
plot filepair('0.1024') u 1:(  1**pow*@LOZ),\
     filepair('0.0512') u 1:(  2**pow*@LOZ),\
     filepair('0.0256') u 1:(  4**pow*@LOZ),\
     filepair('0.0128') u 1:(  8**pow*@LOZ),\
     filepair('0.0064') u 1:( 16**pow*@LOZ),\
     filepair('0.0032') u 1:( 32**pow*@LOZ),\
     filepair('0.0016') u 1:( 64**pow*@LOZ),\
     filepair('0.0008') u 1:(128**pow*@LOZ)


pow=2
set title 'Check at NLO'
set label 1 NLOZ
plot filepair('0.1024') u 1:(  1**pow*@NLOZ),\
     filepair('0.0512') u 1:(  2**pow*@NLOZ),\
     filepair('0.0256') u 1:(  4**pow*@NLOZ),\
     filepair('0.0128') u 1:(  8**pow*@NLOZ),\
     filepair('0.0064') u 1:( 16**pow*@NLOZ),\
     filepair('0.0032') u 1:( 32**pow*@NLOZ),\
     filepair('0.0016') u 1:( 64**pow*@NLOZ),\
     filepair('0.0008') u 1:(128**pow*@NLOZ)


pow=3
set title 'Check at NNLO'
set label 1 NNLOZ
plot filepair('0.1024') u 1:(  1**pow*@NNLOZ),\
     filepair('0.0512') u 1:(  2**pow*@NNLOZ),\
     filepair('0.0256') u 1:(  4**pow*@NNLOZ),\
     filepair('0.0128') u 1:(  8**pow*@NNLOZ),\
     filepair('0.0064') u 1:( 16**pow*@NNLOZ),\
     filepair('0.0032') u 1:( 32**pow*@NNLOZ),\
     filepair('0.0016') u 1:( 64**pow*@NNLOZ),\
     filepair('0.0008') u 1:(128**pow*@NNLOZ)

set title sprintf('Check of artificial bad NNLO (with NNLO piece scaled by %f)',badfactor)
set label 1 NNLObadZ
plot filepair('0.1024') u 1:(  1**pow*@NNLObadZ),\
     filepair('0.0512') u 1:(  2**pow*@NNLObadZ),\
     filepair('0.0256') u 1:(  4**pow*@NNLObadZ),\
     filepair('0.0128') u 1:(  8**pow*@NNLObadZ),\
     filepair('0.0064') u 1:( 16**pow*@NNLObadZ),\
     filepair('0.0032') u 1:( 32**pow*@NNLObadZ),\
     filepair('0.0016') u 1:( 64**pow*@NNLObadZ),\
     filepair('0.0008') u 1:(128**pow*@NNLObadZ)
}
set output
}
