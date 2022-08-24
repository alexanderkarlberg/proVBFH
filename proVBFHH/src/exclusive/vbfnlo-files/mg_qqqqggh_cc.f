c These are the charged current processes.
c
      SUBROUTINE SQQ_QQGGH_CC(pbar,fsign,qbar,gsign,ph,ANS)
c  Last modified: 1 Sept 2008
c  Terrance figy
c
c adapted for POWHEG by Barbara Jaeger: Feb. 2013
c
c  ANS(0) = sum_{i=1,NCOLOR} ANS(i)
c     ANS(i) = weight(i) * ANS(0)
c     weight(i) = |A(i)|^2/sum_{j=1,NCOLOR} |A(j)|^2
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u s -> d c g g h  
C  
C Crossing   1 is u s -> d c g g h  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      include "genpsr.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      integer nexternal
      parameter(nexternal = 7)
c     sign factors 
      integer fsign(4),gsign(2)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
c      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
      real*8 ans(0:6)
      real*8 pbar(0:3,4),qbar(0:3,2),ph(0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T(0:6)
c
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I, ICOL
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      character*79         hel_buff
      common/to_helicity/  hel_buff

      integer          isum_hel
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      SAVE yfrac, igood, IDUM, jhel
      DATA NGRAPHS /   24/               
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   3),IHEL=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   5),IHEL=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,   7),IHEL=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,   9),IHEL=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  11),IHEL=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  13),IHEL=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  15),IHEL=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  17),IHEL=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  19),IHEL=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  21),IHEL=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  23),IHEL=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  25),IHEL=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  27),IHEL=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  29),IHEL=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  31),IHEL=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  33),IHEL=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  35),IHEL=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  37),IHEL=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  39),IHEL=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  41),IHEL=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  43),IHEL=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  45),IHEL=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  47),IHEL=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  49),IHEL=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  51),IHEL=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  53),IHEL=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  55),IHEL=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  57),IHEL=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  59),IHEL=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  61),IHEL=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  63),IHEL=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  72/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      IPROC=1
       
      do ICOL=0,6
      ANS(ICOL) = 0D0
      enddo
      
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
         
         DO IHEL=1,NCOMB
            IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
C               write(*,*) 'IHEL',IHEL
               call CC_QQ_QQGGH(pbar,qbar,ph,nhel(1,ihel),fsign,gsign,T)
               do ICOL=0,6
                  ANS(ICOL) = ANS(ICOL) + T(ICOL)
               enddo
               IF (T(0) .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                  GOODHEL(IHEL,IPROC)=.TRUE.
                  NGOOD = NGOOD +1
                  IGOOD(NGOOD) = IHEL
               ENDIF
            ENDIF
         ENDDO
         JHEL = 1
         ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE                      !RANDOM HELICITY
         DO J=1,ISUM_HEL
            JHEL=JHEL+1
            IF (JHEL .GT. NGOOD) JHEL=1
            HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
            IHEL = IGOOD(JHEL)
            call CC_QQ_QQGGH(pbar,qbar,ph,nhel(1,ihel),fsign,gsign,T)
            do ICOL=0,6
               ANS(ICOL)=ANS(ICOL)+T(ICOL)*HWGT
            enddo
         ENDDO
      ENDIF
               
c     ANS = ANS /DBLE(IDEN(IPROC))
c     now color and spin sum is returned
c
      return
      END
       
       
      subroutine CC_QQ_QQGGH(pbar,qbar,ph,nhel,fsign,gsign,matrix)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u s -> d c g g h  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  6) 
      include "genpsr.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  36, NCOLOR=   6) 
      integer nexternal
      parameter(nexternal = 7)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      real*8 pbar(0:3,4),qbar(0:3,2),ph(0:3)
      real*8 matrix(0:NCOLOR),lsum
      INTEGER NHEL(NEXTERNAL)
      integer fsign(4),gsign(2)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(6,NWAVEFUNCS),WX1(6),WX2(6)
c     local switch
      logical lintOff
      common/interference/ lintOff
      real*8 tfac(ncolor)
C  
C GLOBAL VARIABLES
C  
      include "coupl.inc"
      include 'nlegborn.h'
      include 'pwhg_flst.h' 
      include 'pwhg_flst_2.h' 
      include 'tags.h'
      integer tags(nlegreal)
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,6  ) /    16,    0,   -2,    0,    2,    2/    
C               T[3,1,5,6]T[4,2]                                           
      DATA Denom(2  )/            1/                                       
      DATA (CF(i,2  ),i=1  ,6  ) /     0,   16,    0,    2,    0,    0/    
C               T[3,1,5]T[4,2,6]                                           
      DATA Denom(3  )/            1/                                       
      DATA (CF(i,3  ),i=1  ,6  ) /    -2,    0,   16,    0,    2,    2/    
C               T[3,1,6,5]T[4,2]                                           
      DATA Denom(4  )/            1/                                       
      DATA (CF(i,4  ),i=1  ,6  ) /     0,    2,    0,   16,    0,    0/    
C               T[3,1,6]T[4,2,5]                                           
      DATA Denom(5  )/            1/                                       
      DATA (CF(i,5  ),i=1  ,6  ) /     2,    0,    2,    0,   16,   -2/    
C               T[3,1]T[4,2,5,6]                                           
      DATA Denom(6  )/            1/                                       
      DATA (CF(i,6  ),i=1  ,6  ) /     2,    0,    2,    0,   -2,   16/    
C               T[3,1]T[4,2,6,5]                                           
C ----------
C BEGIN CODE
C ----------                        
      call vetographs(fsign,gsign,tfac)
c      tfac=1
      if (realequiv_tag * alr_tag > 0) then 
         stop 'mg_qqqqggh_cc: both positive' 
      endif

      if (realequiv_tag > 0) then 
         tags = flst_realtags(:,realequiv_tag)
      elseif (alr_tag > 0) then 
         tags = flst_alrtags(:,alr_tag)
      endif

C     kill amplitudes where gluons are connected to wrong fermion lines 
C     (based on comments above about colour structure of the amplitudes) 
      if (any(tags(1:2) == 0) .or. any(tags(4:7) == 0)) then 
C     running without fulltags, so do nothing here 
      elseif (flst_realgluontags(1) == 1 .and. 
     C        flst_realgluontags(2) == 1) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         tfac(2) = 0d0 
         tfac(4:6) = 0 
      elseif (flst_realgluontags(1) == 2 .and. 
     C        flst_realgluontags(2) == 2) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         tfac(1:4) = 0d0 
      elseif (flst_realgluontags(1) == 1 .and. 
     C        flst_realgluontags(2) == 2) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         tfac(1) = 0d0 
         tfac(3:6) = 0 
      elseif (flst_realgluontags(1) == 2 .and. 
     C        flst_realgluontags(2) == 1) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         tfac(1:3) = 0d0 
         tfac(5:6) = 0 
      else
         write(*,*) 'tags', tags 
         stop 'mg_qqqqggh_cc: invalid tags'
      endif

c      JAMP(   1) = -AMP(   1)+AMP(   5)-AMP(  13)+AMP(  17)-AMP(  18)
c      JAMP(   2) = -AMP(   2)-AMP(   4)-AMP(  14)-AMP(  16)
c      JAMP(   3) = -AMP(   3)-AMP(   5)-AMP(   6)-AMP(  15)-AMP(  17)
c      JAMP(   4) = -AMP(   7)-AMP(   9)-AMP(  19)-AMP(  21)
c      JAMP(   5) = -AMP(   8)+AMP(  11)-AMP(  20)+AMP(  23)-AMP(  24)
c      JAMP(   6) = -AMP(  10)-AMP(  11)-AMP(  12)-AMP(  22)-AMP(  23)


      CALL IXXXXX(Pbar(0,1),ZERO ,NHEL(1)*fsign(1),fsign(1),W(1,1))        
      CALL IXXXXX(Pbar(0,3),ZERO ,NHEL(3)*fsign(3),fsign(3),W(1,2))        
      CALL OXXXXX(Pbar(0,2),ZERO ,NHEL(2)*fsign(2),fsign(2),W(1,3))        
      CALL OXXXXX(Pbar(0,4),ZERO ,NHEL(4)*fsign(4),fsign(4),W(1,4))        
      CALL VXXXXX(qbar(0,1),ZERO ,NHEL(5),gsign(1),W(1,5))        
      CALL VXXXXX(qbar(0,2),ZERO ,NHEL(6),gsign(2),W(1,6))        
      CALL SXXXXX(Ph,+1,W(1,7))         
c
      CALL JIOXXX(W(1,2   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,1   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL JVSXXX(W(1,8   ),W(1,7   ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     10  ))                                                          
      if(tfac(1).ne.0) then
         CALL FVOXXX(W(1,3   ),W(1,10  ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
         CALL IOVXXX(W(1,9   ),W(1,11  ),W(1,5   ),GG ,AMP(1   ))
      endif
      
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL FVIXXX(W(1,2   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      if(tfac(2).ne.0) then
         CALL JIOXXX(W(1,12  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
         CALL JIOXXX(W(1,13  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
         CALL VVSXXX(W(1,14  ),W(1,15  ),W(1,7   ),GWWH ,AMP(2   ))           
      endif
      
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,16  ))     
      if(tfac(3).ne.0) then
         CALL IOVXXX(W(1,12  ),W(1,16  ),W(1,10  ),GWF ,AMP(3   ))            
      endif
      
      CALL FVOXXX(W(1,4   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,17  ))     
      if(tfac(2).ne.0) then
         CALL JIOXXX(W(1,2   ),W(1,17  ),GWF ,WMASS   ,WWIDTH  ,W(1,18  ))    
         CALL VVSXXX(W(1,14  ),W(1,18  ),W(1,7   ),GWWH ,AMP(4   ))           
      endif
      
      CALL JVVXXX(W(1,6   ),W(1,5   ),G ,ZERO    ,ZERO    ,W(1,19  ))      
      if(tfac(1).ne.0.and.tfac(3).ne.0) then
         CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,19  ),GG ,AMP(5   ))             
      endif
      
      if(tfac(3).ne.0) then
         CALL FVIXXX(W(1,12  ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,20  ))     
         CALL IOVXXX(W(1,20  ),W(1,3   ),W(1,10  ),GWF ,AMP(6   ))            
      endif
      
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,21  ))     
      if(tfac(4).ne.0) then
         CALL JIOXXX(W(1,9   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,22  ))    
         CALL JIOXXX(W(1,21  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
         CALL VVSXXX(W(1,22  ),W(1,23  ),W(1,7   ),GWWH ,AMP(7   ))           
      endif
      
      CALL JIOXXX(W(1,1   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,24  ))    
      CALL JVSXXX(W(1,24  ),W(1,7   ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     25  ))                                                          
      if(tfac(5).ne.0) then
         CALL FVOXXX(W(1,4   ),W(1,25  ),GWF ,ZERO    ,ZERO    ,W(1,26  ))    
         CALL IOVXXX(W(1,13  ),W(1,26  ),W(1,5   ),GG ,AMP(8   ))             
      endif
      
      if(tfac(4).ne.0) then
         CALL JIOXXX(W(1,1   ),W(1,16  ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
         CALL VVSXXX(W(1,27  ),W(1,23  ),W(1,7   ),GWWH ,AMP(9   ))           
      endif
      
      if(tfac(6).ne.0) then
         CALL IOVXXX(W(1,21  ),W(1,17  ),W(1,25  ),GWF ,AMP(10  ))            
      endif
      
      if(tfac(5).ne.0.and.tfac(6).ne.0) then
         CALL IOVXXX(W(1,2   ),W(1,26  ),W(1,19  ),GG ,AMP(11  ))             
      endif
      
      if(tfac(6).ne.0) then
         CALL FVIXXX(W(1,21  ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
         CALL IOVXXX(W(1,28  ),W(1,4   ),W(1,25  ),GWF ,AMP(12  ))            
      endif
      
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,29  ))     
      if(tfac(1).ne.0) then
         CALL IOVXXX(W(1,9   ),W(1,29  ),W(1,10  ),GWF ,AMP(13  ))            
      endif
      
      if(tfac(2).ne.0) then
         CALL JIOXXX(W(1,1   ),W(1,29  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
         CALL VVSXXX(W(1,30  ),W(1,15  ),W(1,7   ),GWWH ,AMP(14  ))           
      endif
      
      if(tfac(3).ne.0) then
         CALL FVIXXX(W(1,1   ),W(1,10  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
         CALL IOVXXX(W(1,31  ),W(1,16  ),W(1,5   ),GG ,AMP(15  ))             
      endif
      
      if(tfac(2).ne.0) then
         CALL VVSXXX(W(1,30  ),W(1,18  ),W(1,7   ),GWWH ,AMP(16  ))           
      endif
      
      if(tfac(1).ne.0.and.tfac(3).ne.0) then
         CALL IOVXXX(W(1,31  ),W(1,3   ),W(1,19  ),GG ,AMP(17  ))             
      endif
      
      if(tfac(1).ne.0) then
         CALL FVOXXX(W(1,29  ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,32  ))     
         CALL IOVXXX(W(1,1   ),W(1,32  ),W(1,10  ),GWF ,AMP(18  ))            
      endif

      CALL FVOXXX(W(1,4   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,33  ))     
      if(tfac(4).ne.0) then
         CALL JIOXXX(W(1,2   ),W(1,33  ),GWF ,WMASS   ,WWIDTH  ,W(1,34  ))    
         CALL VVSXXX(W(1,22  ),W(1,34  ),W(1,7   ),GWWH ,AMP(19  ))           
      endif
      
      if(tfac(5).ne.0) then
         CALL IOVXXX(W(1,13  ),W(1,33  ),W(1,25  ),GWF ,AMP(20  ))            
      endif
      
      if(tfac(4).ne.0) then
         CALL VVSXXX(W(1,27  ),W(1,34  ),W(1,7   ),GWWH ,AMP(21  ))           
      endif
      
      if(tfac(6).ne.0) then
         CALL FVIXXX(W(1,2   ),W(1,25  ),GWF ,ZERO    ,ZERO    ,W(1,35  ))    
         CALL IOVXXX(W(1,35  ),W(1,17  ),W(1,5   ),GG ,AMP(22  ))             
      endif
      
      if(tfac(6).ne.0.and.tfac(5).ne.0) then
         CALL IOVXXX(W(1,35  ),W(1,4   ),W(1,19  ),GG ,AMP(23  ))             
      endif
      
      if(tfac(5).ne.0) then
         CALL FVOXXX(W(1,33  ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,36  ))     
         CALL IOVXXX(W(1,2   ),W(1,36  ),W(1,25  ),GWF ,AMP(24  ))            
      endif


      JAMP(   1) = -AMP(   1)+AMP(   5)-AMP(  13)+AMP(  17)-AMP(  18)
      JAMP(   2) = -AMP(   2)-AMP(   4)-AMP(  14)-AMP(  16)
      JAMP(   3) = -AMP(   3)-AMP(   5)-AMP(   6)-AMP(  15)-AMP(  17)
      JAMP(   4) = -AMP(   7)-AMP(   9)-AMP(  19)-AMP(  21)
      JAMP(   5) = -AMP(   8)+AMP(  11)-AMP(  20)+AMP(  23)-AMP(  24)
      JAMP(   6) = -AMP(  10)-AMP(  11)-AMP(  12)-AMP(  22)-AMP(  23)


      if (realequiv_tag * alr_tag > 0) then 
         stop 'mg_qqqqggh_cc: both positive' 
      endif

      if (realequiv_tag > 0) then 
         tags = flst_realtags(:,realequiv_tag)
      elseif (alr_tag > 0) then 
         tags = flst_alrtags(:,alr_tag)
      endif

C     kill amplitudes where gluons are connected to wrong fermion lines 
C     (based on comments above about colour structure of the amplitudes) 
      if (any(tags(1:2) == 0) .or. any(tags(4:7) == 0)) then 
C     running without fulltags, so do nothing here 
      elseif (flst_realgluontags(1) == 1 .and. 
     C        flst_realgluontags(2) == 1) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         JAMP(2) = 0d0 
         JAMP(4:6) = 0 
      elseif (flst_realgluontags(1) == 2 .and. 
     C        flst_realgluontags(2) == 2) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         JAMP(1:4) = 0d0 
      elseif (flst_realgluontags(1) == 1 .and. 
     C        flst_realgluontags(2) == 2) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         JAMP(1) = 0d0 
         JAMP(3:6) = 0 
      elseif (flst_realgluontags(1) == 2 .and. 
     C        flst_realgluontags(2) == 1) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         JAMP(1:3) = 0d0 
         JAMP(5:6) = 0 
      else
         write(*,*) 'tags', tags 
         stop 'mg_qqqqggh_cc: invalid tags'
      endif

CC     kill amplitudes where gluons are connected to wrong fermion lines 
CC     (based on comments above about colour structure of the amplitudes) 
C      if (any(tags(1:2) == 0) .or. any(tags(4:7) == 0)) then 
CC     running without fulltags, so do nothing here 
C      elseif (sum(tags(1:nlegreal)) == 8) then ! case 1 2 -> 1 2 1 1 or permutations of it 
C         JAMP(2) = 0d0 
C         JAMP(4:6) = 0 
C      elseif (sum(tags(1:nlegreal)) == 10) then ! case 1 2 -> 1 2 2 2 or permutations of it 
C         JAMP(1:4) = 0d0 
C      elseif (sum(tags(1:nlegreal)) == 9 .and. 
C     C        tags(6) == 1 .and. tags(7) == 2) then ! case 1 2 -> 1 2 1 2 
C         JAMP(1) = 0d0 
C         JAMP(3:6) = 0 
C      elseif (sum(tags(1:nlegreal)) == 9 .and. 
C     C        tags(6) == 2 .and. tags(7) == 1) then ! case 1 2 -> 1 2 2 1 
C         JAMP(1:3) = 0d0 
C         JAMP(5:6) = 0 
C      else
C         write(*,*) 'tags', tags 
C         stop 'mg_qqqqggh_cc: invalid tags'
C      endif




c
      if(lintOff) then
         cf(2,4) = 0d0
         cf(4,2) = 0d0
         cf(5,1) = 0d0
         cf(1,5) = 0d0
         cf(6,1) = 0d0
         cf(1,6) = 0d0
         cf(5,3) = 0d0
         cf(3,5) = 0d0
         cf(6,3) = 0d0
         cf(3,6) = 0d0
      endif
c     eliminate graphs with V-> q qbar
c     for gsign(1 or 2) = -1
      call vetographs(fsign,gsign,tfac)
c
      matrix(0) = 0.D0 
      lsum = 0.D0
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)*tfac(J)
          ENDDO
          matrix(0) = matrix(0) + ZTEMP*DCONJG(JAMP(I))/DENOM(I)*tfac(I) 
          lsum       = lsum + (DBLE(JAMP(I))**2+DIMAG(JAMP(I))**2)*tfac(I)
      ENDDO
      
      do I=1,NCOLOR
         matrix(I) = matrix(0) * (DBLE(JAMP(I))**2 + DIMAG(JAMP(I))**2) 
     $        * tfac(I)         !/sum
         if(abs(lsum).gt.0.D0) matrix(I)=matrix(I)/lsum
      enddo
c    
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     this subroutine computes tfac(ncolor) 
c     tfac is used to eliminate gauge graphs with V-> q qbar.
c
      subroutine vetographs(fsign,gsign,tfac)
      implicit none 
      include 'nlegborn.h'
      include 'pwhg_flst.h' 
      include 'pwhg_flst_2.h' 
      integer fsign(4),gsign(2),i
      real*8 tfac(6)
      integer tags(nlegreal) 
      logical, save :: first_time = .true. 

      if (first_time) then 
         write(*,*) 'changed vetographs for gg channel' 
         first_time = .false. 
      endif
      if (realequiv_tag > 0) then 
         tags = flst_realtags(:,realequiv_tag)
      elseif (alr_tag > 0) then 
         tags = flst_alrtags(:,alr_tag)
      endif
c
c   by default tfac(i) = one    
      do i=1,6
         tfac(i) = 1d0
      enddo
c
      if((gsign(1).eq.-1).and.(gsign(2).eq.1)) then !one gluon
         if((fsign(3).eq.fsign(4)).and.
     $        (fsign(1).eq.-fsign(2))) then !initial gluon on 21 line
            tfac(5) = 0
            tfac(6) = 0
            tfac(4) = 0
         elseif((fsign(3).eq.-fsign(4)).and.
     $           (fsign(1).eq.fsign(2))) then !initial gluon on 43 line
            tfac(1) = 0
            tfac(2) = 0
            tfac(3) = 0
         endif
         
      elseif((gsign(1).eq.1).and.(gsign(2).eq.-1)) then !one gluon
         if((fsign(3).eq.fsign(4)).and.
     $        (fsign(1).eq.-fsign(2))) then !initial gluon on 21 line    
             tfac(5) = 0
             tfac(6) = 0
             tfac(2) = 0
          elseif((fsign(3).eq.-fsign(4)).and.
     $            (fsign(1).eq.fsign(2))) then !initial gluon on 43 line
             tfac(1) = 0
             tfac(4) = 0
             tfac(3) = 0
          endif
      elseif((gsign(1).eq.-1).and.(gsign(2).eq.-1)) then !two gluons
         if((fsign(3).eq.-fsign(4)).and. !gluon qbar(mu,2) on 43 line
     $        (fsign(1).eq.-fsign(2))) then !gluon qbar(mu,1) on 21 line
            tfac(1) = 0
            tfac(3) = 0
            tfac(5) = 0
            tfac(6) = 0
            tfac(4) = 0

         else
            print*,'error in vetograph'
            print*, gsign(1:2)
            print*, fsign(1:4)
            stop
         endif
      endif
      return
      end
