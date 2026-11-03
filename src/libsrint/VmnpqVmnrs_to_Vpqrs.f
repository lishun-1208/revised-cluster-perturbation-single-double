      SUBROUTINE VmnpqVmnrs_to_Vpqrs(Vmnpq, Vmnrs, SYMV1, SYMV2, 
     &                POPM, POPN, iopmn,POPP, POPQ,ioppq, POPR, POPS,
     &                ioprs, Vpqrs, SYMV3, FACT, fact2)
C This program is to calculate 
C   V(pq,rs)=fact2*V(pq,rs) + fact*V1(mn,pq)*V2(mn,rs). 
C SYMV1, SYMV2 and SYMV3 are the symmetry of Vmnpq, Vmnrs and Vpqrs.
C One should not that
C the order of index in V1 and V2  must be like (mn,pq) and (mn,rs). this
C program can be employed to deal with IIJKL, IIJKA and IABIJ. but for
C IABCI and IABCD, the cost of memory would be very large so that some 
C special treatment must be employed.
C---if iopmn = 0, m<n,  1 m,n   2 m<=n    
C---if ioppq = 0, p<q,  1 p,q   2 p<=q
C---if ioprs = 0, r<s,  1 r,s   2 r<=s
C     
      use mod_symm
      use mod_iop 
      IMPLICIT NONE 
C
      INTEGER, intent(in) :: SYMV1, SYMV2, SYMV3
      INTEGER, intent(in) :: IOPPQ, IOPMN, IOPRS 
      INTEGER, intent(in) :: POPP(8),POPQ(8),POPR(8)
      INTEGER, intent(in) :: POPS(8),POPM(8),POPN(8)
      real*4, intent(in),managed :: Vmnpq(*), Vmnrs(*)
      real*4, intent(in),managed :: Vpqrs(*)
      DOUBLE PRECISION, intent(in) :: FACT, fact2
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
      INTEGER, EXTERNAL :: IDSYMSZ,IRPDSO 
      DOUBLE PRECISION, EXTERNAL :: DNRM2
C
C compute size of pq and rm for each irrep 
C
      if(dabs(fact2)>1.d-10.and.dabs(fact2-1.d0)>1.d-10) then 
         write(6,*)'attention in VV_to_V!, fact2/=0.d0 and fact2/=1.d0'
      endif
      if(dirprd(symv1, symv2)/=symv3) then 
         write(6,*)'error of symmetry VmnpqVmnrs', symv1, symv2, symv3
         stop
      endif
      CALL IZERO(INUMPQ,8); CALL IZERO(INUMRS,8)
      CALL IZERO(INUMMN,8)
C-----inumpq
      
        DO IRREP=1, NIRREP
           INUMPQ(IRREP)  = IRPDSO(IRREP, POPP, POPQ, IOPPQ) 
        ENDDO
C-----inumrs
        DO IRREP=1, NIRREP
              INUMRS(IRREP)  = IRPDSO(IRREP, POPR, POPS, IOPRS) 
        ENDDO
C-----inummn
        DO IRREP=1, NIRREP
              INUMMN(IRREP)  = IRPDSO(IRREP, POPM, POPN, IOPMN) 
        ENDDO
C compute absolute offsets to beginning of each irrep 
C---for Vmnpq
      IABSOFF1(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV1)
         IABSOFF1(IRREPR+1) = IABSOFF1(IRREPR)
     &                 +INUMMN(IRREPL)*INUMPQ(IRREPR)*ISD
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
         IOFFV1 = IABSOFF1(IRREPPQ)
         IOFFV2 = IABSOFF2(IRREPRS) 
         IOFFV3 = IABSOFF3(IRREPRS)
         CALL XGEMMSO('T','N',INUMPQ(IRREPPQ), INUMRS(IRREPRS),  
     &        INUMMN(IRREPMN), FACT, Vmnpq(IOFFV1),INUMMN(IRREPMN),
     &        Vmnrs(IOFFV2),INUMMN(IRREPMN), Fact2,
     &        Vpqrs(IOFFV3),INUMPQ(IRREPPQ))
      ENDDO
C
      RETURN
      END  



















 
