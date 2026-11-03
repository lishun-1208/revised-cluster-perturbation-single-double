      SUBROUTINE VpqrmVpqrn_to_Vmn(Vpqrm, Vpqrn, SYMV1, SYMV2, 
     &                POPP, POPQ,ioppq,  POPR, POPM, POPN,
     &               Vmn, SYMV3, FACT)
C This program is to calculate
C Vmn=Vmn + fact*V2(pq,rm)*V1(pq,rn). 
C One should not that
C the order of index in V2 and V1  must be like (pq,rm) and (pq,rn). 
C r<m or r<n can not be included in the calculation.
C program can be employed to deal with IIJKL, IIJKA and IABIJ. but for
C IABCI and IABCD, the cost of memory would be very large so that some 
C special treatment must be employed.
C
C---if ioppq = 0, p<q,  else if ioppq = 1, p,q =2 p<=q
C
      use mod_iop
      use mod_symm
      IMPLICIT NONE 
C
      INTEGER, intent(in) :: SYMV1, SYMV2, SYMV3
      INTEGER, intent(in) :: IOPpq
      INTEGER, intent(in) :: POPP(8),POPQ(8),POPR(8)
      INTEGER, intent(in) :: POPM(8),POPN(8)
      REAL*4, intent(in) :: Vpqrm(*), Vpqrn(*)
      REAL*4, intent(in) :: Vmn(*)
      DOUBLE PRECISION, intent(in) :: FACT
C
      INTEGER INUMPQ, INUMRM, INUMRN
      INTEGER IABSOFF1(9), IPOSV1(8,8), IABSOFF2(9)
      INTEGER IPOSV2(8,8), IPOSV3(8)
      INTEGER IRR, IRR1, IRR2
      INTEGER IRREPR, IRREPL, IRREP, IRREP1
      INTEGER IRREPPQ, IRREPRS, IRREPMN 
      INTEGER IRREPP, IRREPQ, IRREPR1, IRREPS1, IRREPM, IRREPN
      INTEGER NUML, NUMR, IOFFV1, IOFFV2, IOFFV3
      INTEGER NUMM, NUMN, NUMPQ, NUMRM, NUMRN, IRREPRM, IRREPRN 
      INTEGER IOFFCHECK2(8), IOFFCHECK3(8), IRREP2R 
      INTEGER NSIZE, nump, numq, numr1, nums
      INTEGER, EXTERNAL :: IDSYMSZ, IRPDSO 
      DOUBLE PRECISION, EXTERNAL :: DNRM2
C
C compute size of pq and rm for each irrep 
C
C-----inumpq
C---for Vpqrm
      IABSOFF1(1) = 1
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR, SYMV1)
         NUMR=IRPDSO(IRREPR, POPR, POPM, 1)
         NUML=IRPDSO(IRREPL, POPP, POPQ, IOPPQ)
         IABSOFF1(IRREPR+1) = IABSOFF1(IRREPR)+NUML*NUMR*ISD
         IPOSV1(1, IRREPR) = 1 
         DO IRREP1 = 1, NIRREP-1 
            IRREP2R=DIRPRD(IRREP1, IRREPR)
            IPOSV1(IRREP1+1, IRREPR)=IPOSV1(IRREP1,IRREPR)+
     &              POPR(IRREP2R)*POPM(IRREP1)*ISD
         ENDDO
      ENDDO
C---for Vpqrn
      IABSOFF2(1) = 1
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR, SYMV2)
         NUMR=IRPDSO(IRREPR, POPR, POPN, 1)
         NUML=IRPDSO(IRREPL, POPP, POPQ, IOPPQ)
         IABSOFF2(IRREPR+1) = IABSOFF2(IRREPR)+NUML*NUMR*ISD
         IPOSV2(1, IRREPR) = 1 
         DO IRREP1 = 1, NIRREP-1 
            IRREP2R=DIRPRD(IRREP1, IRREPR)
            IPOSV2(IRREP1+1, IRREPR)=IPOSV2(IRREP1,IRREPR)+
     &              POPR(IRREP2R)*POPN(IRREP1)*ISD
         ENDDO
      ENDDO
C---for Vmn
      IPOSV3(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV3)
         NUMR = POPN(IRREPR)
         NUML = POPM(IRREPL)
         IPOSV3(IRREPR+1)=IPOSV3(IRREPR)+NUMR*NUML*ISD
      ENDDO
C
      DO IRREPRM = 1, NIRREP
         IRREPPQ = DIRPRD(IRREPRM, SYMV1)
         IRREPRN = DIRPRD(IRREPPQ, SYMV2)
         NUMRM = IRPDSO(IRREPRM, POPR, POPM, 1) 
         NUMPQ = IRPDSO(IRREPPQ, POPP, POPQ, IOPPQ) 
         DO IRREPM = 1, NIRREP 
            IRREPR = DIRPRD(IRREPM, IRREPRM)
            IRREPN = DIRPRD(IRREPRN, IRREPR) 
            IRREPMN = DIRPRD(IRREPM, IRREPN) 
            IF(IRREPMN/=SYMV3) THEN 
               WRITE(6,*)'ERROR IN SYMMETRY'
               write(6,*)'symmry', irrepmn, symv1, symv2, symv3
               stop
            endif
            NUMM=POPM(IRREPM)
            NUMR=POPR(IRREPR)
            NUMN=POPN(IRREPN)   
           IOFFV1 = IABSOFF1(IRREPRM)+
     &            (IPOSV1(IRREPM, IRREPRM)-1)*NUMPQ
           IOFFV2 = IABSOFF2(IRREPRN)+
     &            (IPOSV2(IRREPN, IRREPRN)-1)*NUMPQ
           IOFFV3 = IPOSV3(IRREPN)
           CALL XGEMMSO('T','N',numm,numn,numpq*numr,  
     &               FACT, Vpqrm(IOFFV1),numpq*numr,
     &               Vpqrn(IOFFV2),numpq*numr, 1.d0,
     &               Vmn(IOFFV3),numm)
         ENDDO
      ENDDO
C
      RETURN
      END  



















 
