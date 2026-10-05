      program asq
      implicit none
      real*8 m, a, alphasPDF
      external alphasPDF
      call initPDFsetByName('NNPDF30_nnlo_as_0118')
      call initPDF(0)
 10   read(*,*,end=20) m, a
      write(*,'(a,f9.3,a,f9.5,a,f9.5,a,f8.4)') ' mur', m, '  powheg', a,
     $     '  lhapdf', alphasPDF(m), '  ratio', a/alphasPDF(m)
      goto 10
 20   end
