      SUBROUTINE VpqrmVms_to_Vpqrs(Vpqrm, Vms, SYMV2, SYMV1, 
     &                POPM, POPP, POPQ, IOP, POPR, POPS,
     &                Vpqrs, SYMV3, FACT, FACT2)
C This program is to calculate V2(p,q,r,m)*V1(m,s). One should not that
C the order of index in V2 and V1  must be like (p,q,r,m) and (m,s). this
C program can be employed to deal with IIJKL, IIJKA and IABIJ. but for
C IABCI and IABCD, the cost of memory would be very large so that some 
C special treatment must be employed.
C      
C----the choice of IOP (i,j,k,l) means occupied while a,b,c,d means
C    virtual orbitals)
C 
C----IOP = 0 means P<Q,R,S 
C----IOP = 1 means P,Q,R,S
C
C----Vpqrs = FACT2*Vpqrs + FACT*Vpqrm*Vms
C
      use mod_iop
      use mod_symm
      IMPLICIT NONE 
C
      INTEGER, intent(in) :: SYMV1, SYMV2, SYMV3, IOP
      INTEGER, intent(in) :: POPP(8),POPQ(8),POPR(8),POPS(8),POPM(8)  
      real*4, intent(in),managed :: Vpqrm(*), Vms(*)
      real*4, intent(in),managed :: Vpqrs(*)
      real*8, intent(in) :: FACT, FACT2
C
      INTEGER INUMLFT(8), INUMRHT(8), INUMRHT3(8), INUMRHT32(8)
      INTEGER IABSOFF(9), IPOSLFT(8,8), IPOSRHT(8,8), IPOSV1(8)
      INTEGER IABSOFF3(9), IPOSRHT3(8,8), IABSOFF32(9)
      INTEGER IRREP, IRR1,IRR2, IRREPL, IRREPR, IRREPR3
      INTEGER IRREP1, IRREP2L, IRREP2R, IRREPRS
      INTEGER IRREPM, IRREPQ, IRREPP,  IRREPS1, IRREPR1
      INTEGER NUML, NUMR, IOFFV1, IOFFV2, IOFFV3
      INTEGER IOFFCHECK2(8), IOFFCHECK3(8) 
      INTEGER NSIZE, nump, numq, numr1, nums, numm
      INTEGER, EXTERNAL :: IDSYMSZ, IRPDSO 
      DOUBLE PRECISION, EXTERNAL :: DNRM2
C
C compute size of pq and rm for each irrep 
C
      if(dabs(fact2)>1.d-10.and.dabs(fact2-1.d0)>1.d-10) then 
         write(6,*)'attention in VV_to_V!, fact2/=0.d0 and fact2/=1.d0'
      endif
      if(dirprd(symv1, symv2)/=symv3) then 
         write(6,*)'error of symmetry VpqrmVms', symv1, symv2, symv3
         stop
      endif
      IF(IOP/=0.AND.IOP/=1.AND.IOP/=2) THEN 
         WRITE(6,*)'ERROR OF IOP IN VpqrmVsm'
         stop
      endif
C-----inum=p,q,r,s
      CALL IZERO(INUMLFT,8); CALL IZERO(INUMRHT,8)
      CALL IZERO(INUMRHT3,8); CALL IZERO(INUMRHT32,8)
      DO IRREP=1, NIRREP
         INUMLFT(IRREP) = IRPDSO(IRREP, POPP, POPQ, IOP)
         INUMRHT(IRREP) = IRPDSO(IRREP, POPR, POPM, 1)
         INUMRHT3(IRREP) = IRPDSO(IRREP, POPR, POPS, 1)
      ENDDO
C
C compute absolute offsets to beginning of each irrep 
C---for Vpqrm
      CALL IZERO(IABSOFF,8)
      IABSOFF(1) = 1
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR, SYMV2)
         IABSOFF(IRREPR+1) = IABSOFF(IRREPR)+INUMLFT(IRREPL)*INUMRHT(IRREPR)*ISD
         IPOSRHT(1,IRREPR)=1
         DO IRREP1 = 1, NIRREP-1
            IRREP2R=DIRPRD(IRREP1,IRREPR)
            IPOSRHT(IRREP1+1,IRREPR)=IPOSRHT(IRREP1,IRREPR)+
     &                            POPR(IRREP2R)*POPM(IRREP1)*ISD
         ENDDO             
      ENDDO
C---for Vpqrs
      CALL IZERO(IABSOFF3,8)
      CALL IZERO(IABSOFF32,8)
      IABSOFF3(1) = 1; IABSOFF32(1) = 1
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR, SYMV3)
         IABSOFF3(IRREPR+1) = 
     &        IABSOFF3(IRREPR)+INUMLFT(IRREPL)*INUMRHT3(IRREPR)*ISD
         IABSOFF32(IRREPR+1) =
     &        IABSOFF32(IRREPR)+INUMLFT(IRREPL)*INUMRHT32(IRREPR)*ISD
         IPOSRHT3(1,IRREPR)=1
         DO IRREP1 = 1, NIRREP-1
            IRREP2R=DIRPRD(IRREP1,IRREPR)
            IRREP2L=DIRPRD(IRREP1,IRREPL)
            IPOSRHT3(IRREP1+1,IRREPR)=IPOSRHT3(IRREP1,IRREPR)+
     &                            POPR(IRREP2R)*POPS(IRREP1)*ISD
         ENDDO
      ENDDO
      
C---for Vms
      IPOSV1(1) = 1 
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV1) 
         NUMR = POPS(IRREPR)
         NUML = POPM(IRREPL)
         IPOSV1(IRREPR+1)=IPOSV1(IRREPR)+NUMR*NUML*ISD
      ENDDO
C---check the number 
C
      DO IRREPR = 1, NIRREP 
         IRREPL = DIRPRD(IRREPR, SYMV2) 
         IRREPRS = DIRPRD(IRREPL, SYMV3)
         DO IRREPS1 = 1, NIRREP 
            IRREPM = DIRPRD(IRREPS1, SYMV1)
            IRREPR1 = DIRPRD(IRREPM, IRREPR) 
            IOFFV1 = IPOSV1(IRREPS1)
            IOFFV2 = IABSOFF(IRREPR)+
     &           (IPOSRHT(IRREPM, IRREPR)-1)*INUMLFT(IRREPL)
            IOFFV3 = IABSOFF3(IRREPRS)+
     &           (IPOSRHT3(IRREPS1, IRREPRS)-1)*INUMLFT(IRREPL)
            CALL XGEMMSO('N','N',INUMLFT(IRREPL)*POPR(IRREPR1), 
     &        POPS(IRREPS1), POPM(IRREPM), FACT, 
     &        Vpqrm(IOFFV2),INUMLFT(IRREPL)*POPR(IRREPR1),
     &        Vms(IOFFV1),POPM(IRREPM), FACT2,
     &        Vpqrs(IOFFV3),INUMLFT(IRREPL)*POPR(IRREPR1))
         ENDDO
      ENDDO
C
      RETURN
      END  



















 
