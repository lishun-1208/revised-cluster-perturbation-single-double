      SUBROUTINE VmnpqVmnrs_to_Vpqrs_large_2(Vmnpq, Vmnrs, SYMV1, SYMV2,
     &                POPM, POPN,iopmn, POPP, POPQ,ioppq, POPR,POPS,
     &                ioprs, Vpqrs, SYMV3, FACT,fact2, Vmnrsc, icore)
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
      real*4, intent(inout) :: icore(*)
      real*4, intent(inout) :: Vmnrsc(*)
      real*4, device, intent(inout) :: Vmnpq(*)
      real*4, device, intent(inout) :: Vpqrs(*), Vmnrs(*)
      real*4, intent(in) :: FACT, fact2
      double precision :: tt1, tt2, rt1, rt2, copytime
C      character(len=8), intent(in) :: filename1
C
      real*4 test
      INTEGER INUMPQ(8), INUMRS(8), INUMMN(8)
      INTEGER IABSOFF1(8), IABSOFF3(8), IABSOFF2(8)
      INTEGER IRR, IRR1, IRR2, iopv1, iopv2, iopv3
      INTEGER IRREPR, IRREPL, IRREP
      INTEGER IRREPPQ, IRREPRS, IRREPMN 
      INTEGER IRREPP, IRREPQ, IRREPR1, IRREPS1, IRREPM, IRREPN
      INTEGER*8 NUML, NUMR, IOFFV1, IOFFV2, IOFFV3
      INTEGER IOFFCHECK2(8), IOFFCHECK3(8) 
      INTEGER NSIZE, nump, numq, numr1, nums, numm
      INTEGER numn, n, iopf, nummn, numpq, numrs, iter
      INTEGER nfirst, v3start, zorb, nread, nreadold 
      INTEGER*8 MEMSIZE, nread8
      INTEGER*8 USIZE, i, j
      integer  nleft, nsizev2, nsizev3, istat
      INTEGER, EXTERNAL :: IDSYMSOC, irpdso, isymoffso 
C
      zorb = 1
      USIZE = 1024
      MEMSIZE = USIZE*USIZE*USIZE/4 !4GB
C      write(6,*)'memsize', memsize
      nsizev2 = idsymsoc(symv2, popm, popn, iopmn, popr, pops, ioprs)
      nsizev3 = idsymsoc(symv3, popp, popq, ioppq, popr, pops, ioprs)
C
C compute size of pq and rm for each irrep 
C
       iopf = 1 
C      if(ioppq==0.and.iopmn==0) iopf=1
C      if(ioppq==1.and.iopmn==1) iopf=2
C      if(ioppq==2.and.iopmn==2) iopf=1
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
     &                      +INUMMN(IRREPL)*INUMPQ(IRREPR)
      ENDDO
C---for Vmnrs
      IABSOFF2(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV2)
         IABSOFF2(IRREPR+1) = IABSOFF2(IRREPR)
     &                      +INUMMN(IRREPL)*INUMRS(IRREPR)
      ENDDO
C---for Vpqrs
      IABSOFF3(1) = 1
      DO IRREPR = 1, NIRREP-1
         IRREPL = DIRPRD(IRREPR, SYMV3)
         IABSOFF3(IRREPR+1) = IABSOFF3(IRREPR)
     &                      +INUMPQ(IRREPL)*INUMRS(IRREPR)
      ENDDO
C
      copytime = 0.d0
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
C         write(6,*)'nread', memsize, nread
C         stop
         iter = 1
         DO
            write(6,*)'iter', iter 
            IF(NREAD>NUMRS) NREAD=NUMRS
            if(nleft==0) exit 
!           if(iter==1) then 
!               call myicopy_gr(vmnrsc(ioffv2), vmnrs, nread*nummn)
!           else 
!               call cublasscopy(nread*nummn, vmnrs(memsize+10), 1, vmnrs, 1)
!           endif
                nread8=nread
                ioffv2 = ioffv2 + nread8*nummn
                nleft =  nleft - nread 
                nfirst = nfirst+ nread
                nreadold = nread
                if(nread>nleft) nread = nleft
             if(nleft>0) then 
C!$omp parallel 
C!$omp sections
C!$omp section
!               call cublassgemm('T','N', numpq, nreadold,
!    &           nummn, fact, Vmnpq, nummn, Vmnrs, nummn,
!    &           fact2, Vpqrs(ioffv3),inumpq(irreppq))
C!$omp section
                call myicopy_gr(vmnrsc(ioffv2),vmnrs(memsize+10), nread*nummn)
C!$omp end sections
C!$omp end parallel 
             elseif(nleft==0) then 
!               call cublassgemm('T','N', numpq, nreadold,
!    &           nummn, fact, Vmnpq, nummn, Vmnrs, nummn,
!    &           fact2, Vpqrs(ioffv3),inumpq(irreppq))
             endif
                ioffv3 = ioffv3 + nreadold*numpq
                iter = iter+1
!$omp parallel 
!$omp sections
!$omp section
C                do i = 1, 20000
C                  write(6,*)'i', i
C                enddo

!$omp section 
C                do j = 1, 20000
C                  write(6,*)'j', j
C                enddo

!$omp end sections
!$omp end parallel 
         ENDDO     
      ENDDO
C
      RETURN
      END  



















 
