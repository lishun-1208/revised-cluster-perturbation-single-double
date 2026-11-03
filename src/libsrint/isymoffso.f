
      INTEGER FUNCTION ISYMOFFSO(IRREPQ1,IRREP, POPP, POPQ, IOP) 
C
C THIS ROUTINE RETURNS THE start point of SIZE OF V(PQ) with special
C IRREPQ
C    iop==0, p<q
C    iop==1, p,q
C    iop==2, p<=q
C
CEND
      use mod_symm
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION POPP(8), POPQ(8)
      integer, external :: irpdso
C
      ISYMOFFSO = 1
      if(iop==1) then  
      DO 10 IRREPQ=1,irrepq1-1
       IRREPP=DIRPRD(IRREPQ,IRREP)
       ISYMOFFSO=ISYMOFFSO+
     & POPP(IRREPP)*POPQ(IRREPQ)
10    CONTINUE
      elseif(iop==0) then
C       write(6,*)'test in isymoffso' 
       do irrepq = 1, irrepq1-1
          irrepp = dirprd(irrepq, irrep) 
C       write(6,*)'test in isymoffso2',  
C     &         popp(irrepp), popp(irrepq)
          if(irrepp<irrepq) then 
             isymoffso = isymoffso + 
     &           POPP(IRREPP)*POPQ(IRREPQ)
C       write(6,*)'test in isymoffso3', isymoffso 
          elseif(irrepp==irrepq) then 
             isymoffso = isymoffso + 
     &           POPP(IRREPP)*(POPP(IRREPQ)-1)/2
C       write(6,*)'test in isymoffso4', isymoffso 
          endif
        enddo
      elseif(iop==2) then 
       do irrepq = 1, irrepq1-1
          irrepp = dirprd(irrepq, irrep) 
          if(irrepp<irrepq) then 
             isymoffso = isymoffso + 
     &           POPP(IRREPP)*POPQ(IRREPQ)
          elseif(irrepp==irrepq) then 
             isymoffso = isymoffso + 
     &           POPP(IRREPP)*(POPP(IRREPQ)+1)/2
          endif
        enddo
      endif
C
      RETURN
      END
