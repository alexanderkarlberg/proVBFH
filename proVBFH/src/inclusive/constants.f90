module constants
  use helper
  integer , parameter :: maxloops = 5
  real(dp), parameter :: pi = 3.1415926535897932384626433_dp, pi2 = pi**2
  real(dp), parameter :: pi3 = pi**3, pi4 = pi**4, pi5 = pi**5, pi6 = pi**6
  real(dp), parameter :: zeta2 = pi2 / 6.0_dp, zeta3 = 1.2020569031595942853997_dp
  real(dp), parameter :: zeta4 = pi4 / 90.0_dp, zeta5 = 1.03692775514336992633136_dp
  real(dp), parameter :: zeta6 = pi6 / 945.0_dp, zeta7 = 1.00834927738192282683979_dp
  real(dp), parameter :: EulerGamma = 0.57721566490153286061_dp
  ! These numbers are needed for 4-loop decoupling. Can be found in hep-ph/0512060
  real(dp), parameter :: T54_3 = -8445.8046390310298_dp, T62_2 = -4553.4004372195263_dp
  real(dp), parameter :: B4=-1.7628000870737708640618976346798188072151372743890_dp
  real(dp), parameter :: A4=0.517479061673899386330758161898862945622377475141_dp
  real(dp), parameter :: A5=0.50840057924226870745910884925858994131954112566_dp

end module constants
