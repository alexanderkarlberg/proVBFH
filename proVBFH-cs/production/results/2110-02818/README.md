# Stable-Higgs VBF fiducial cross sections of arXiv:2110.02818 with proVBFH-cs

Card: `../../2110-02818/{powheg.input,vbfnlo.input}` (13 TeV, NNPDF31_nnlo_as_0118 = LHAPDF id 303600,
mu_R = mu_F = mu0(pT,H) with factors 1/2, 1, 2 together, jets anti-kt R = 0.4, pT >= 25 GeV, |y| <= 4.5,
tagging jets = the two leading jets of that list, |dy| > 4.5, mjj > 600 GeV, opposite hemispheres,
M_H = 125, Gamma_H = 4.165 MeV (narrow-width Higgs, so it does not enter), M_W = 80.398, Gamma_W = 2.105,
M_Z = 91.1876, Gamma_Z = 2.4952, G_F = 1.16639e-5, scheme (G_F, M_W, M_Z), diagonal CKM).
Binary: proVBFH-cs-p1506 built from this tree (copy in /ptmp/mpp/akarlber/cs2110/bin). Results are
plain means over seeds (`summarize.py`: error of the mean from the seed scatter), >= 2-jet histogram
`sig(all VBF cuts 2 jets)`, in fb. `<order>-<J|H>-W<k>.top` are all histograms of the combined runs
(W1 mu0, W2 mu0/2, W3 2mu0); J = selection of 2110.02818, H = selection of 1506.02660.
`numbers.txt` is the output of `summarize.py`.

## Set-up changes with respect to the p1506 card
- PDF: `lhans 261000` is NNPDF30_nnlo_as_0118 in LHAPDF's index; 2110.02818 needs NNPDF31 (303600; set
  installed in ~/.local/share/LHAPDF). alpha_s is taken from the PDF at all orders.
- G_F: the card had 1.16637e-5, and the inclusive part hard-coded it (incl_parameters.f90:11). It now reads
  FERMI_CONST (default unchanged); the card has 1.16639e-5.
- W width 2.141 -> 2.105, H width 4.030e-3 -> 4.165e-3 (card only; the Higgs width is unused at `higgsbreitwigner 0`).
- Analysis: new flag `jet_hardest` (p1506_analysis.f; 0 default = jets restricted to pT, |y| first, then the
  two leading ones = 2110.02818; 1 = the two hardest jets of the event must pass pT and |y| = 1506.02660).
  At LO the two give identical results (two partons).

## Comparison (fb); the paper quotes integers, so allow about 0.5 fb for rounding
| order | scale | proVBFH-cs | paper | diff | pull |
|---|---|---|---|---|---|
| LO (100 jobs) | mu0 | 970.58 +- 0.45 | 971 | -0.4 | -0.9 |
| | mu0/2 | 1039.57 +- 0.49 | 1040 | -0.4 | -0.9 |
| | 2mu0 | 909.18 +- 0.42 | 910 | -0.8 | -2.0 (-0.8 with rounding) |
| NLO (40 jobs) | mu0 | 890.09 +- 0.68 | 890 | +0.1 | +0.1 |
| | mu0/2 | 872.50 +- 0.69 | 872 | +0.5 | +0.7 |
| | 2mu0 | 897.30 +- 0.67 | 898 | -0.7 | -1.0 |
| NNLO (100 jobs, binary before the virtual fix) | mu0 | 854.10 +- 1.07 | 859 | -4.9 | -4.6 |
| | mu0/2 | 843.45 +- 1.32 | 849 | -5.5 | -4.2 |
| | 2mu0 | 864.20 +- 0.91 | 867 | -2.8 | -3.1 |

LO and NLO agree at all three scales. NNLO is 0.3-0.6% below the paper (this is the deficit later traced to the
H+3j virtual, see below).

## Selection effect (1506-style minus 2110-style is negative; J - H shown, same seeds)
| order | mu0 | mu0/2 | 2mu0 |
|---|---|---|---|
| LO | 0 | 0 | 0 |
| NLO (40 seeds) | +0.887 +- 0.008 | +1.166 +- 0.010 | +0.694 +- 0.007 |
| NNLO (40 seeds) | +1.20 +- 0.16 | +1.28 +- 0.21 | +1.09 +- 0.12 |

Taking the two hardest jets overall (1506) instead of the two leading jets within acceptance lowers the cross
section by 0.1% at NLO and 0.14% at NNLO. It does not explain the NNLO deficit.

## Cost
LO 100 jobs about 0.1 core-h, NLO 80 jobs about 16 core-h, NNLO-J 100 jobs (+12 reruns after node failures on gt37 and gt28)
about 236 core-h, NNLO-H 40 jobs about 85 core-h.
