      subroutine gen_uborn_idx
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      integer j,alr,em
      do j=1,rad_nkinreg
         rad_kinreg_on(j)=.false.
      enddo
      call pick_random(flst_nborn,rad_btilde_arr,rad_ubornidx)
      rad_alr_nlist=flst_born2alr(0,rad_ubornidx)
      do j=1,rad_alr_nlist
         alr=flst_born2alr(j,rad_ubornidx)
         em=flst_emitter(alr)
         rad_alr_list(j)=alr
         if(em.le.2) then
            rad_kinreg_on(1)=.true.
         else
            rad_kinreg_on(em-flst_lightpart+2)=.true.
         endif
      enddo
      end

      subroutine gen_real_idx
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      call pick_random(rad_alr_nlist,rad_real_arr,rad_realidx)
      rad_realalr=rad_alr_list(rad_realidx)
      kn_emitter = flst_emitter(rad_realalr)
      end

      subroutine pick_random(n,r,jret)
      implicit none
      integer n,jret
      real * 8 r(n)
      real * 8 tot(0:n)
      integer j
      real * 8 ran,random
      tot(0)=0
      do j=1,n
         tot(j)=tot(j-1)+r(j)
      enddo
      ran=tot(n)*random()
      do j=1,n
         if(tot(j).gt.ran) then
            jret=j
            return
         endif
      enddo
      write(*,*) ' ***********************************'
      write(*,*) ' POWHEGBOX:pick_random: could not pick a value'
      write(*,*) ' set choice to 1'
      write(*,*) ' ***********************************'
      jret=1
c      call exit(-1)
      end
      

      subroutine gen_remnant(iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      integer iret
      real * 8 dum1,dum2,dum3,dum4
      real * 8 random
      external random
c Remnant or regular?
      if(random().gt.rad_reg_tot/(rad_reg_tot+rad_damp_rem_tot)) then
c remnant
         iret=1
         call pick_random(flst_nalr,rad_damp_rem_arr,rad_realalr)
         rad_ubornidx=flst_alr2born(rad_realalr)
         kn_emitter=flst_emitter(rad_realalr)
         if(kn_emitter.le.2) then
            call gen_real_phsp_isr(rad_xradremn,dum1,dum2,dum3,dum4)
         else
            call gen_real_phsp_fsr(rad_xradremn,dum1,dum2,dum3)
         endif
      else
c regular
         iret=2
         call pick_random(flst_nregular,rad_reg_arr,rad_realreg)
         call gen_real_phsp_isr(rad_xradremn,dum1,dum2,dum3,dum4)
      endif
      end
