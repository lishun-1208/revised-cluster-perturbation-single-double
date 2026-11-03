      SUBROUTINE iassymso(Vin, Vout, SYMV,
     &                POPP, POPQ, POPR, POPS, IOP)
C
C----This program is to exchange Vpqrs----
C----iop=0, Vout(p<q,r<s) = Vin(p<q,r,s)-Vin(p<q,s,r)
C----iop=1, Vout(p,q,r<s) = Vin(p,q,r,s)-Vin(p,q,s,r)
C----iop=2, Vout(p,q,r,s) = Vin(p,q,r,s)-Vin(q,p,s,r)
C----iop=3, Vout(p,q,r,s) = Vin(p,q,r,s)+Vin(q,p,s,r)
C----iop=4, Vout(p<q,r<s) = Vin(p,q,r<s)-Vin(q,p,r<s)
C----iop=5, Vout(p<q,r,s) = Vin(p,q,r,s)-Vin(q,p,r,s)
C 
      use mod_iop
      use mod_symm
      IMPLICIT NONE 
C
      INTEGER, intent(in) :: SYMV, IOP
      INTEGER, intent(in) :: POPP(8),POPQ(8),POPR(8),POPS(8)  
      real*4, intent(in) :: Vin(*)
      real*4, intent(in) :: Vout(*)
C
      INTEGER INUMLFT(8), INUMRHT(8), INUMRHT2(8), INUMLFT2(8)
      INTEGER IABSOFF(8), IABSOFF4(8) 
      INTEGER IABSOFF3(8), IPOSRHT3(8,8), IABSOFF2(8)
      INTEGER IRREP, IRR1,IRR2, IRREPL, IRREPR, IRREPR3
      INTEGER IRREP1, IRREP2L, IRREP2R, IRREPRS
      INTEGER IRREPM, IRREPQ, IRREPP,  IRREPS1, IRREPR1
      INTEGER NUML, NUMR, IOFFV1, IOFFV2, IOFFV3
      INTEGER IOFFCHECK2(8), IOFFCHECK3(8) 
      INTEGER NSIZE, nump, numq, numr1, nums, numm
      INTEGER ISTAT
      INTEGER, EXTERNAL :: IDSYMSZ, IRPDSO 
      DOUBLE PRECISION, ALLOCATABLE :: SCR(:,:)
      DOUBLE PRECISION, EXTERNAL :: DNRM2
C
C compute size of pq and rm for each irrep 
C
C-----inum=p,q,r,s
      CALL IZERO(INUMLFT,8); CALL IZERO(INUMRHT,8)
      CALL IZERO(INUMLFT2,8); CALL IZERO(INUMRHT2,8)
      DO IRREP=1, NIRREP
         INUMLFT(IRREP) = IRPDSO(IRREP, POPP, POPQ, 1)
         INUMRHT(IRREP) = IRPDSO(IRREP, POPR, POPS, 1)
         INUMLFT2(IRREP) = IRPDSO(IRREP, POPP, POPQ, 0)
         INUMRHT2(IRREP) = IRPDSO(IRREP, POPR, POPS, 0)
      ENDDO
C
C---for Vpqrs
      CALL IZERO(IABSOFF,8)
      CALL IZERO(IABSOFF2,8)
      IABSOFF(1) = 1; IABSOFF2(1) = 1
      IABSOFF3(1) = 1; IABSOFF4(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV)
C---start of V(p,q,r,s)
         IABSOFF(IRREPR+1) = 
     &        IABSOFF(IRREPR)+INUMLFT(IRREPL)*INUMRHT(IRREPR)*ISD
C---start of V(p,q,r<s)
         IABSOFF2(IRREPR+1) =
     &        IABSOFF2(IRREPR)+INUMLFT(IRREPL)*INUMRHT2(IRREPR)*ISD
C---start of V(p<q,r,s)
         IABSOFF3(IRREPR+1) = 
     &        IABSOFF3(IRREPR)+INUMLFT2(IRREPL)*INUMRHT(IRREPR)*ISD
C---start of V(p<q,r<s)
         IABSOFF4(IRREPR+1) = 
     &        IABSOFF4(IRREPR)+INUMLFT2(IRREPL)*INUMRHT2(IRREPR)*ISD
      ENDDO
C    
      if(iop==1) then 
C----iop=0, Vout(p,q,r<s) = Vin(p,q,r,s)-Vin(p,q,s,r)
        do irrepr = 1, nirrep 
           irrepl = dirprd(irrepr, symv) 
           call assym2so(irrepr, popr, inumlft(irrepl),
     &                Vout(iabsoff2(irrepr)), Vin(iabsoff(irrepr))) 
        enddo
      elseif(iop==0) then
C----iop=1, Vout(p<q,r<s) = Vin(p<q,r,s)-Vin(p<q,s,r)
        do irrepr = 1, nirrep 
           irrepl = dirprd(irrepr, symv) 
           call assym2so(irrepr, popr, inumlft2(irrepl),
     &                Vout(iabsoff4(irrepr)), Vin(iabsoff3(irrepr))) 
        enddo
C
      elseif(iop==2) then
C----iop=2, Vout(p,q,r,s) = Vin(p,q,r,s)-Vin(q,p,s,r) 
        do irrepr = 1, nirrep
           irrepl = dirprd(irrepr, symv)
           numl = inumlft(irrepl)
           call dasymrhfso(irrepl, irrepr, popp, popr, numl,
     &               Vout(iabsoff(irrepr)), Vin(iabsoff(irrepr)))
        enddo
C
      elseif(iop==3) then 
C----iop=3, Vout(p,q,r,s) = Vin(p,q,r,s)+Vin(q,p,s,r) 
        do irrepr = 1, nirrep 
           irrepl = dirprd(irrepr, symv)
           numl = inumlft(irrepl)
           call dsymrhfso(irrepl, irrepr, popp, popr, numl,
     &               Vout(iabsoff(irrepr)), Vin(iabsoff(irrepr)))
        enddo
      elseif(iop==4) then 
C----iop=4, Vout(p<q,r<s) = Vin(p,q,r<s)-Vin(q,p,r<s)
        do irrepr = 1, nirrep 
           irrepl = dirprd(irrepr, symv) 
           call assym3so(irrepl, popp, inumlft2(irrepl),
     &                inumlft(irrepl), inumrht2(irrepr),
     &                Vout(iabsoff4(irrepr)), Vin(iabsoff2(irrepr)))
        enddo
      elseif(iop==5) then
C----iop=4, Vout(p<q,r,s) = Vin(p,q,r,s)-Vin(q,p,r,s)
        do irrepr = 1, nirrep
           irrepl = dirprd(irrepr, symv)
           call assym3so(irrepl, popp, inumlft2(irrepl),
     &                inumlft(irrepl), inumrht(irrepr),
     &                Vout(iabsoff3(irrepr)), Vin(iabsoff(irrepr)))   
        enddo
      else
         write(6,*)'error of iop in iassymso'
         stop
      endif 
C
      RETURN
      END  



















 
