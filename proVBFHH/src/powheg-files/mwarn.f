      subroutine mwarn(str)
      character *(*) str
      write(*,*) '********** WARNING **********'
      write(*,*) '*********  ',str,'  *********'
      end
