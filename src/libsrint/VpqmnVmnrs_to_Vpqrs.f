      SUBROUTINE VpqmnVmnrs_to_Vpqrs(Vpqmn, Vmnrs, SYMV1, SYMV2, 
     &                POPM, POPN,iopmn, POPP, POPQ,ioppq, POPR,POPS,
     &                ioprs, Vpqrs, SYMV3, FACT,fact2)
C This program is to calculate
C Vpqrs=fact2*Vpqrs + fact*V2(pq,mn)*V1(mn,rs). 
C One should not that
C the order of index in V2 and V1  must be like (pq,mn) and (mn,rs). this
C program can be employed to deal with IIJKL, IIJKA and IABIJ. but for
C IABCI and IABCD, the cost of memory would be very large so that some 
C special treatment must be employed.
C
C---if iopmn = 0, m<n, 1 m,n 2 m<=n      
C---if ioppq = 0, p<q, 1 p,q 2 p<=q
C---if ioprs = 0, r<s, 1 r,s 2 r<=s
C
      use mod_iop
      use mod_symm
      IMPLICIT NONE 
C
      INTEGER, intent(in) :: SYMV1, SYMV2, SYMV3
      INTEGER, intent(in) :: IOPpq, iopmn, ioprs
      INTEGER, intent(in) :: POPP(8),POPQ(8),POPR(8)
      INTEGER, intent(in) :: POPS(8),POPM(8),POPN(8)
      REAL*4, intent(in),managed :: Vpqmn(*), Vmnrs(*)
      REAL*4, intent(in),managed :: Vpqrs(*)
      REAL*8, intent(in) :: FACT, fact2
C
      INTEGER INUMPQ(8), INUMRS(8), INUMMN(8)
      INTEGER IABSOFF1(8), IABSOFF3(8), IABSOFF2(8)
      INTEGER IRR, IRR1, IRR2
      INTEGER IRREPR, IRREPL, IRREP
      INTEGER IRREPPQ, IRREPRS, IRREPMN 
      INTEGER IRREPP, IRREPQ, IRREPR1, IRREPS1, IRREPM, IRREPN
      INTEGER NUML, NUMR, IOFFV1, IOFFV2, IOFFV3
      INTEGER IOFFCHECK2(8), IOFFCHECK3(8) 
      INTEGER NSIZE, nump, numq, numr1, nums, numm
      INTEGER, EXTERNAL :: IDSYMSZ, irpdso 
      DOUBLE PRECISION, EXTERNAL :: DNRM2
C
C compute size of pq and rm for each irrep 
C
      if(dabs(fact2)>1.d-10.and.dabs(fact2-1.d0)>1.d-10) then 
         write(6,*)'attention in VV_to_V!, fact2/=0.d0 and fact2/=1.d0'
      endif
      if(dirprd(symv1, symv2)/=symv3) then 
         write(6,*)'error of symmetry in vpqmnvmnrs', symv1, symv2, symv3
         stop
      endif
      CALL IZERO(INUMPQ,8); CALL IZERO(INUMRS,8)
      CALL IZERO(INUMMN,8)
C-----inumpq
       do irrep=1, nirrep
          inumpq(irrep) = irpdso(irrep, popp, popq, ioppq) 
       enddo
C-----inumrs
       do irrep=1, nirrep
          inumrs(irrep) = irpdso(irrep, popr, pops, ioprs) 
       enddo
C-----inummn
       do irrep=1, nirrep
          inummn(irrep) = irpdso(irrep, popm, popn, iopmn) 
       enddo
C compute absolute offsets to beginning of each irrep 
C---for Vpqmn
      IABSOFF1(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV1)
         IABSOFF1(IRREPR+1) = IABSOFF1(IRREPR)
     &                 +INUMPQ(IRREPL)*INUMMN(IRREPR)*ISD
      ENDDO
C---for Vmnrs
      IABSOFF2(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV2)
         IABSOFF2(IRREPR+1) = IABSOFF2(IRREPR)
     &                 +INUMMN(IRREPL)*INUMRS(IRREPR)*ISD
      ENDDO
C---for Vpqrs
      IABSOFF3(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV3)
         IABSOFF3(IRREPR+1) = IABSOFF3(IRREPR)
     &                 +INUMPQ(IRREPL)*INUMRS(IRREPR)*ISD
      ENDDO
      
      DO IRREPRS = 1, NIRREP 
         IRREPPQ = DIRPRD(IRREPRS, SYMV3) 
         IRREPMN = DIRPRD(IRREPPQ, SYMV1)
         IOFFV1 = IABSOFF1(IRREPMN)
         IOFFV2 = IABSOFF2(IRREPRS) 
         IOFFV3 = IABSOFF3(IRREPRS)
         CALL XGEMMSO('N','N',INUMPQ(IRREPPQ), INUMRS(IRREPRS),  
     &        INUMMN(IRREPMN), FACT, Vpqmn(IOFFV1),INUMPQ(IRREPPQ),
     &        Vmnrs(IOFFV2),INUMMN(IRREPMN), fact2,
     &        Vpqrs(IOFFV3),INUMPQ(IRREPPQ))
      ENDDO
C
      RETURN
      END  



















 
