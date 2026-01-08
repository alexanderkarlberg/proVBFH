#!/usr/bin/env python
from math import *
asOld=0.112639920858 # at 125 GeV
asNew=0.125156764869 # at 62.5 GeV

# results at xmuR=xmuF=2
LO=0.888190981382; NLO=0.919476466783; NNLO=0.914825262835

nloCoeff = (NLO-LO)/asOld
nnloCoeff = (NNLO-NLO)/asOld**2

print nloCoeff, nnloCoeff
# gives 0.277747757302 -0.36659008089

NLOnew = LO + nloCoeff*asNew
print NLOnew
# gives 0.922952992136

b0=0.61
nnloCoeffNew = nnloCoeff - b0*log(4.0)*nloCoeff
NNLOnew = NLOnew + nnloCoeffNew*asNew**2
print NNLOnew
# gives 0.913531521284
