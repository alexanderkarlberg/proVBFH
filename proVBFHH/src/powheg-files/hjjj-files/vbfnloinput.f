      function vbfnloinput(stringa)
      implicit none
      real * 8 vbfnloinput
      character *(*) stringa
      integer maxnum
      parameter (maxnum=150)
      character * 100 line,line0
      character * 20 string
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 20 keywords(maxnum)
      real * 8 values(maxnum)
      logical used(maxnum)
      integer ios,numvalues,j,k,l,imode
      integer ini
      data ini/0/      
      save ini, keywords, values, numvalues,used
      string=stringa
      if(ini.eq.0) then
         open(unit=33,file='vbfnlo.input',status='old',iostat=ios)
         if(ios.ne.0) then
            write(*,*) ' Enter the prefix for this run < 20 characters'
            read(*,'(a)') pwgprefix
            do lprefix=20,1,-1
               if(pwgprefix(lprefix:lprefix).ne.' ') then
                  goto 11
               endif
            enddo
 11         continue
            lprefix=lprefix+1
            if(lprefix.gt.20) lprefix=20
            pwgprefix(lprefix:lprefix)='-'
            open(unit=33,file=pwgprefix(1:lprefix)//'vbfnlo.input',
     #           status='old',iostat=ios)
            if(ios.ne.0) then            
               write(*,*) ' cannot open ',
     #              pwgprefix(1:lprefix)//'vbfnlo.input'
               stop
            endif
         else
            pwgprefix='pwg'
            lprefix=3
         endif
         numvalues=0
         do l=1,1000000
            line0=' '
            read(unit=33,fmt='(a)',iostat=ios) line0
            if(ios.ne.0.and.line0.eq.' ') goto 10
            line=line0
            do k=1,100
               if(line(k:k).eq.'#'.or.line(k:k).eq.'!') then
                  line(k:)=' '
               endif
            enddo
            if(line.ne.' ') then
               if(numvalues.eq.maxnum) then
                  write(*,*) ' too many entries in vbfnloinput.dat'
                  call exit(-1)
               endif
               numvalues=numvalues+1
c skip blanks
 12            if(line(1:1).eq.' ') then
                  line=line(2:)
                  goto 12
               endif
               k=index(line,' ')
               keywords(numvalues)=line(1:k-1)
               line=line(k+1:)
               read(unit=line,fmt=*,iostat=ios) values(numvalues)
               used(numvalues)=.false.
               if(ios.ne.0) then
                  write(*,*) ' vbfnloinput error: cannot parse '
                  write(*,'(a)') line0
                  stop
               endif
            endif
         enddo
 10      continue
         close(33)
         ini=1
      endif
      if(stringa.eq.'print unused tokens') then
         do j=1,numvalues
            if(.not.used(j)) then
               write(*,*)'vbfnloinput WARNING: unused variable ',
     1              keywords(j)
            endif
         enddo
         return
      endif
      if(string(1:1).eq.'#') then
         string=string(2:)
         imode=0
      else
         imode=1
      endif
      do j=1,numvalues
         if(string.eq.keywords(j)) then
            vbfnloinput=values(j)
            if(.not.used(j)) then
               used(j)=.true.
               write(*,*) ' vbfnloinput keyword ',keywords(j),
     1                    ' set to ',values(j)
            endif
            return
         endif
      enddo
      if(imode.eq.1) then
         write(*,*) ' vbfnloinput: keyword ',string,' not found'
         call exit(-1)
      endif
c Not found; assign value -1d6; store the token anyhow
      if(numvalues.eq.maxnum) then
         write(*,*) ' too many entries in vbfnloinput.dat'
         write(*,*) ' increase maxnum in vbfnloinput.f'
         call exit(-1)
      endif
      numvalues=numvalues+1
      keywords(numvalues)=string
      values(numvalues)=-1d6
      used(numvalues)=.true.
      vbfnloinput=-1d6
      write(*,*) ' vbfnloinput keyword ',keywords(j),
     1     ' absent; set to ',values(j)
      end
