      program pdft
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      real*8 pdf(-22:22), f(-6:6), x, q
      integer i
      call init_flsttag
      call init_phys
      x = 0.05d0
      q = 100d0
      st_mufact2 = q**2
      call pdfcall(1, x, pdf)
      call initPDFsetByName('NNPDF30_nnlo_as_0118')
      call initPDF(0)
      call evolvePDF(x, q, f)
      do i = -3, 3
         write(*,'(a,i3,a,es14.6,a,es14.6,a,f10.5)') ' f', i, ' powheg ',
     $        pdf(i), ' lhapdf f(x) ', f(i)/x, ' ratio ', pdf(i)*x/f(i)
      enddo
      write(*,*) 'st_alpha', st_alpha, ' st_lambda5MSB', st_lambda5MSB
      end
