*****************************************************
***                                               ***
***  Declaration of model parameters.             ***
***  sophy@particle.uni-karlsruhe.de              ***
***                                               ***
*****************************************************

*      double precision Pi
*      parameter (Pi = 3.1415926535897932384626433832795029D0)

** complex i
	double complex cI
	parameter (cI = (0D0,1D0))


** em coupling constant, electrix charge, fermi constant
      double precision AlfaQED, EL, gfermi, AlfaQED2, delAlfa


** Electroweak mixing angle and fermion and guage boson masses and their squares
** as well as fermion mass array
      double precision SW, SW2, CW, CW2
      double precision MW, MZ, MW2, MZ2

      double precision MsU, MD, MS, MC, MB, MT
      double precision MU2, MD2, MS2, MC2, MB2, MT2

      double precision ME, MM, ML
      double precision ME2, MM2, ML2

      double precision Mf(4,3), Mf2(4,3)

** Higgs mixing angles, masses, propagator matrix and width
      double precision MH, MH2

** Coupling factors: the values that the charges and left/right axial vector 
** couplings can take
      double precision Qu, Qd, gL3, gL4, gR3, gR4, dgL3, dgL4, dgR3, dgR4


** Renormalisation parameters: mass^2 of photon, mass of quark, UV divergence 
      double precision massG, massQ, Divergence


      common /alfaparamsD/ AlfaQED, EL, AlfaQED2, gfermi, delAlfa

      common /swparamsD/ SW, SW2, CW, CW2
      common /mwzparamsD/ MW, MZ, MW2, MZ2

      common /mqparamsD/ MsU, MD, MS, MC, MB, MT
      common /mq2paramsD/ MU2, MD2, MS2, MC2, MB2, MT2
      common /mlparamsD/ ME, MM, ML, ME2, MM2, ML2
      common /mfparamsD/ Mf, Mf2
