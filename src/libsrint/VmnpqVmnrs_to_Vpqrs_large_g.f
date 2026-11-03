      SUBROUTINE VmnpqVmnrs_to_Vpqrs_large_g(Vmnpq, Vmnrs, SYMV1, SYMV2,
     &                POPM, POPN,iopmn, POPP, POPQ,ioppq, POPR,POPS,
     &                ioprs, Vpqrs, SYMV3, FACT,fact2, Vmnrsc)
C This program is to calculate
C Vpqrs=fact2*Vpqrs + fact*V1(mnpq)*V2(mn,rs). 
C with V2(mnrs) be very large like IABCD
C
C---if iopmn = 0, m<n,  else m,n      
C---if ioppq = 0, p<q,  else p,q
C---if ioprs = 0, r<s,  else r,s
C 
      use mod_iop
      use mod_symm
      IMPLICIT NONE 
C
      INTEGER, intent(in) :: SYMV1, SYMV2, SYMV3
      INTEGER, intent(in) :: IOPpq, iopmn, ioprs
      INTEGER, intent(in) :: POPP(8),POPQ(8),POPR(8)
      INTEGER, intent(in) :: POPS(8),POPM(8),POPN(8)
      real*4, intent(inout),device :: Vmnpq(*)
      real*4, intent(inout),device :: Vpqrs(*), Vmnrs(*)
      real*4, intent(in) :: Vmnrsc(*)
      real*8, intent(in) :: FACT, fact2
      double precision :: tt1, tt2, rt1, rt2, copytime
C
      real*4 test
      INTEGER INUMPQ(8), INUMRS(8), INUMMN(8)
      INTEGER IABSOFF1(8), IABSOFF3(8), IABSOFF2(8)
      INTEGER IRR, IRR1, IRR2, iopv1, iopv2, iopv3
      INTEGER IRREPR, IRREPL, IRREP
      INTEGER IRREPPQ, IRREPRS, IRREPMN 
      INTEGER IRREPP, IRREPQ, IRREPR1, IRREPS1, IRREPM, IRREPN
      INTEGER NUML, NUMR, IOFFV1, IOFFV2, IOFFV3
      INTEGER IOFFCHECK2(8), IOFFCHECK3(8) 
      INTEGER NSIZE, nump, numq, numr1, nums, numm
      INTEGER numn, n, iopf, nummn, numpq, numrs
      INTEGER nfirst, v3start, zorb, nread 
      INTEGER MEMSIZE, nread8
      INTEGER USIZE
      integer  nleft, nsizev2, nsizev3, istat
      INTEGER, EXTERNAL :: IDSYMSOC, irpdso, isymoffso 
C
      zorb = 1
      USIZE = 1024
      MEMSIZE = USIZE*USIZE*USIZE*2/4 !2GB
      nsizev2 = idsymsoc(symv2, popm, popn, iopmn, popr, pops, ioprs)
      nsizev3 = idsymsoc(symv3, popp, popq, ioppq, popr, pops, ioprs)
C
C compute size of pq and rm for each irrep 
C
      iopv2 = 1; iopv3 = 1
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
C---for Vmnpq
      IABSOFF1(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV1)
         IABSOFF1(IRREPR+1) = IABSOFF1(IRREPR)
     &                      +INUMMN(IRREPL)*INUMPQ(IRREPR)*ISD
      ENDDO
C---for Vmnrs
      IABSOFF2(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV2)
         IABSOFF2(IRREPR+1) = IABSOFF2(IRREPR)
     &                      +INUMMN(IRREPL)*INUMRS(IRREPR)*ISD
      ENDDO
C---for Vpqrs
      IABSOFF3(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV3)
         IABSOFF3(IRREPR+1) = IABSOFF3(IRREPR)
     &                      +INUMPQ(IRREPL)*INUMRS(IRREPR)*ISD
      ENDDO
C
      DO IRREPRS = 1, NIRREP 
         IRREPPQ = DIRPRD(IRREPRS, SYMV3) 
         IRREPMN = DIRPRD(IRREPPQ, SYMV1)
         NUMPQ = IRPDSO(IRREPPQ, POPP, POPQ, IOPPQ) 
         NUMRS = IRPDSO(IRREPRS, POPR, POPS, IOPRS) 
         NUMMN = IRPDSO(IRREPMN, POPM, POPN, IOPMN)
         IOFFV1 = IABSOFF1(IRREPMN)
         IOFFV2 = IABSOFF2(IRREPRS) 
         IOFFV3 = IABSOFF3(IRREPRS)
         NFIRST = 1
         NLEFT = NUMRS
         NREAD = MEMSIZE/NUMMN
         DO 
            IF(NREAD>NUMRS) NREAD=NUMRS
            if(nleft==0) exit
                call myicopy_g(vmnrsc(ioffv2), vmnrs, nread*nummn)
                call xgemmso('T','N', numpq, nread,
     &           nummn, fact, Vmnpq(ioffv1), nummn, Vmnrs, nummn,
     &           fact2, Vpqrs(ioffv3),numpq)
                ioffv3 = ioffv3 + nread*numpq*isd
                ioffv2 = ioffv2 + nread*nummn*isd
                nleft =  nleft - nread 
                nfirst = nfirst+ nread
                if(nread>nleft) nread = nleft
         ENDDO     
      ENDDO
C
      RETURN
      END  



















 
