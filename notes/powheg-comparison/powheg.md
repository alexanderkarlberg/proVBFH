# How to get agreement with POWHEG

To get agreement with powheg one needs to modify a few files (here we give ther recipe for the changes in powheg)

## init_coupling.f

Powheg uses a fixed alphaem and sinthetaW. Around line 56 paste the following:

```
   ph_Zmass  = 91.1976d0
   ph_Zwidth =  2.4952d0
   ph_Wmass  = 80.377d0
   ph_Wwidth =  2.141d0
   
   ph_sthw2 = 1d0 - (ph_Wmass/ph_Zmass)**2
   ph_alphaem = 1.16638d-5 * sqrt(2d0) * ph_Zmass**2 * ph_sthw2 *
        $     (1d0 - ph_sthw2) / pi
```

Remebering to also update the values for the H,W,Z masses and widths.

## Born_phsp.f

Here we need to set the `BW = .false.` to force narrow width. (line 23)

For a fixed scale run on needs to set `dynamical = .false.` and make sure that the scales are set to

```
   muf=ph_Hmass
   mur=ph_Hmass
```

## Included here

Here we also include the input files from proVBFH and powheg that
gives agreement. At the inclusive level the following command line
should reproduce the numbers are

```
   ../provbfh_incl -lo -sqrts 13600 -pdf NNPDF40MC_nlo_as_01180 -scale-choice 0 -ncall1 1000000 -ncall2 10000000 -zwidth 0.0 -mz 91.1976 -mw 80.377 -wwidth 0.0 -hwidth 0.00407
```

Gives an LO cross section of 4.554737E+00 +- 7.304425E-04