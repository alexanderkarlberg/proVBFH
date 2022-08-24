      subroutine calc_als
      implicit none
c
c define values for the factorization scales and alphas in vbfnlo format
c
      include 'pwhg_st.h'
      include 'scales.inc'
c
      integer L
      real*8 qsq
c
cccccccccccccccccccccccccccccccccccccccccccccc
c
      qsq = st_muren2
      do L = 1,3
         als(1,L) =  st_alpha
         als(2,L) =  st_alpha       
         mursq(1,L) = qsq
	 mursq(2,L) = qsq
      enddo         

      end
