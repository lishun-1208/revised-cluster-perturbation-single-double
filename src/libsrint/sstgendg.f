       SUBROUTINE SSTGENDG(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
C
C THIS ROUTINE RESORTS A W(PQ,RS) SYMMETRY-PACKED QUANTITY TO
C W(PS,RQ) [INTERCHANGE OF Q AND S].  
C
C INPUT :
C         WIN - THE W(PQ,RS) QUANTITY [ALL IRREPS]
C       NSIZE - THE SIZE OF THE W QUANTITY
C        POP1 - POPULATION VECTOR FOR P
C        POP2 - POPULATION VECTOR FOR Q
C        POP3 - POPULATION VECTOR FOR R
C        POP4 - POPULATION VECTOR FOR S
C      IRREPX - THE OVERALL SYMMETRY OF THE LIST (USUALLY 1)
C
C OUTPUT:
C        WOUT - THE RESORTED P(PS,RQ) QUANTITY
C
CEND
      use mod_symm
      IMPLICIT NONE 
C
      CHARACTER*4, intent(in) :: TYPE
      DOUBLE PRECISION, intent(in),device  :: WIN(*)
      DOUBLE PRECISION, intent(out),device :: WOUT(*)
      DOUBLE PRECISION, intent(out) :: IPOS(*)
      INTEGER, intent(in) :: POP1(8),POP2(8),POP3(8),POP4(8), IRREPX
      INTEGER*8 INUMLFTN(9),INUMRHTN(9),IABSOFFN(9)
      INTEGER*8 IPOSRHTN(9,9),IPOSLFTN(9,9), NSIZE
      INTEGER*8 INUMLFTO(9),INUMRHTO(9),IABSOFFO(9)
      INTEGER*8 IPOSRHTO(9,9),IPOSLFTO(9,9), IOLD, INEW
      INTEGER*8 ISTARTN(9,9,9,9), istarto(9,9,9,9)
      INTEGER*8 NLFTO, NLFTN 
C
      integer irrep12, irrep34, irrep14, irrep23, irrep13, irrep24
      integer irrep42, irrep4, irrep3, irrep2, irrep1, irrep41
      integer num4, num3, num2, num1
      integer i4, i3, i2, i1, irrep21, ilfto, ilftn
      integer irrrht, irrlft, irreprht, irreplft, irrep31
      integer irr1, irr2, irrep2r, irrep2l, irrep32, irrep
      integer*8 ioffold, ioffnew
      integer*8 ioffo1, ioffo2, ioffo3, ioffo4, ioffoo
      integer*8 ioffn1, ioffn2, ioffn3, ioffn4, ioffnn
C
      INTEGER, EXTERNAL :: IDSYMSOC
      DOUBLE PRECISION :: DNRM2
C
      IF(TYPE.EQ.'1423')THEN
C
C            P Q R S -> P S Q R
C            1 2 3 4    1 4 2 3 
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8) 
       call izero(inumrhtn,8)
       call izero(inumlfto,8) 
       call izero(inumrhto,8)
       do irrep = 1, nirrep 
          do irr1 = 1, nirrep 
             irr2 = dirprd(irr1, irrep) 
             inumlftn(irrep) = inumlftn(irrep)+pop1(irr2)*pop4(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop2(irr2)*pop3(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
       enddo 
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0 
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep 
          irrlft = dirprd(irrrht, irrepx)  
          iabsoffn(irrrht+1) = iabsoffn(irrrht) + 
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) + 
     &                    inumlfto(irrlft)*inumrhto(irrrht) 
           iposrhto(1, irrrht) = 0 
           iposlfto(1, irrlft) = 0 
           iposrhtn(1, irrrht) = 0 
           iposlftn(1, irrlft) = 0 
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep 
              irrep2r = dirprd(irrep1, irrrht) 
              irrep2l = dirprd(irrep1, irrlft) 
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) + 
     &                                     pop2(irrep2r)*pop3(irrep1)*
     &                                     inumlftn(irrlft) 
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) + 
     &                                     pop1(irrep2L)*pop4(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) + 
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft) 
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) + 
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1, irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep 
          irreplft = dirprd(irreprht, irrepx) 
          do irrep3 = 1, nirrep 
             irrep2 = dirprd(irreprht, irrep3) 
             irrep23 = dirprd(irrep2, irrep3) 
             num3 = pop3(irrep3) 
             num2 = pop2(irrep2) 
                do irrep4 =1, nirrep 
                   irrep1 =dirprd(irreplft, irrep4)
                   irrep14 = dirprd(irrep1, irrep4)
                   irrep12 = dirprd(irrep1, irrep2) 
                   irrep34 = dirprd(irrep3, irrep4) 
                   num4 = pop4(irrep4) 
                   num1 = pop1(irrep1)  
                 iold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + iposlfto(irrep2, irrep12)
                 nlfto = inumlfto(irrep12)
                 inew = iabsoffn(irrep23)+iposrhtn(irrep3, irrep23)
     &                 + iposlftn(irrep4, irrep14)
                 nlftn = inumlftn(irrep14)
!$acc kernels
!$acc loop independent
                   do i3 = 1, num3  
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i4 = 1, num4 
!$acc loop independent
                   do i1 = 1, num1
                 wout(inew+
     &                  ((i3-1)*num2+i2-1)*nlftn
     &                 +  (i4-1)*num1 + i1)
     &            =win(iold+
     &                  ((i4-1)*num3+i3-1)*nlfto
     &                 +  (i2-1)*num1 + i1) 
C                ioffold = iold+ 
C    &                  ((i4-1)*num3+i3-1)*nlfto
C    &                 +  (i2-1)*num1 + i1
C
C                ioffnew =inew+
C    &                  ((i3-1)*num2+i2-1)*nlftn
C    &                 +  (i4-1)*num1 + i1
C                wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo 
C     
      ELSEIF(TYPE.EQ.'1432')THEN
C
C            P Q R S -> P S R Q
C            1 2 3 4    1 4 3 2
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop1(irr2)*pop4(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop3(irr2)*pop2(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep) 
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop2(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop4(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
C
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep14 = dirprd(irrep1, irrep4)
                   irrep32 = dirprd(irrep3, irrep2)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
C             write(6,*)'i', iabsoffn(1), iposrhtn(1,1), inumlftn(1),
C     &                   iposlftn(1,1)
C                 write(6,*)'num4 and num1', num1, num4 
C                 write(6,*)'ipos', iposlftn(irrep4, irrep14), 
C     &                                irrep4, irrep14
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep32)+iposrhtn(irrep2, irrep32)
     &                 + ((i2-1)*num3+i3-1)*inumlftn(irrep14)
     &                 + iposlftn(irrep4, irrep14) + (i4-1)*num1 + i1
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
C
      ELSEIF(TYPE.EQ.'2314')THEN
C
C            P Q R S -> Q R P S
C            1 2 3 4    2 3 1 4
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop2(irr2)*pop3(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop1(irr2)*pop4(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop1(irrep2r)*pop4(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop2(irrep2L)*pop3(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep23 = dirprd(irrep2, irrep3)
                   irrep14 = dirprd(irrep1, irrep4)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
                 iold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + iposlfto(irrep2, irrep12) 
C
                 inew = iabsoffn(irrep14)+iposrhtn(irrep4, irrep14)
     &                 + iposlftn(irrep3, irrep23)
                 nlfto = inumlfto(irrep12)
                 nlftn = inumlftn(irrep23)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 wout(inew
     &                 + ((i4-1)*num1+i1-1)*nlftn
     &                 + (i3-1)*num2 + i2) = 
     &           win(iold
     &                 + ((i4-1)*num3+i3-1)*nlfto
     &                 + (i2-1)*num1 + i1)
C                ioffold = iold
C    &                 + ((i4-1)*num3+i3-1)*nlfto
C    &                 + (i2-1)*num1 + i1
C
C                ioffnew = inew
C    &                 + ((i4-1)*num1+i1-1)*nlftn
C    &                 + (i3-1)*num2 + i2
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
C                wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
      ELSEIF(TYPE.EQ.'4123')THEN
C
C            P Q R S -> S P Q R
C            1 2 3 4    4 1 2 3 
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop4(irr2)*pop1(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop2(irr2)*pop3(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop2(irrep2r)*pop3(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop4(irrep2L)*pop1(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep41 = dirprd(irrep4, irrep1)
                   irrep23 = dirprd(irrep2, irrep3)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep23)+iposrhtn(irrep3, irrep23)
     &                 + ((i3-1)*num2+i2-1)*inumlftn(irrep41)
     &                 + iposlftn(irrep1, irrep41) + (i1-1)*num4 + i4
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'2431')THEN
C
C            P Q R S -> Q S R P
C            1 2 3 4    2 4 3 1 
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop2(irr2)*pop4(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop3(irr2)*pop1(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop1(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop2(irrep2L)*pop4(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep24 = dirprd(irrep2, irrep4)
                   irrep31 = dirprd(irrep1, irrep3)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep31)+iposrhtn(irrep1, irrep31)
     &                 + ((i1-1)*num3+i3-1)*inumlftn(irrep24)
     &                 + iposlftn(irrep4, irrep24) + (i4-1)*num2 + i2
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'4231')THEN
C
C            P Q R S -> S Q R P
C            1 2 3 4    4 2 3 1 
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop4(irr2)*pop2(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop3(irr2)*pop1(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop1(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop4(irrep2L)*pop2(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep42 = dirprd(irrep2, irrep4)
                   irrep31 = dirprd(irrep1, irrep3)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep31)+iposrhtn(irrep1, irrep31)
     &                 + ((i1-1)*num3+i3-1)*inumlftn(irrep42)
     &                 + iposlftn(irrep2, irrep42) + (i2-1)*num4 + i4
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'4132')THEN
C
C            P Q R S -> S P R Q
C            1 2 3 4    4 1 3 2
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop4(irr2)*pop1(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop3(irr2)*pop2(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep) 
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop2(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop4(irrep2L)*pop1(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
C
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep41 = dirprd(irrep1, irrep4)
                   irrep32 = dirprd(irrep3, irrep2)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep32)+iposrhtn(irrep2, irrep32)
     &                 + ((i2-1)*num3+i3-1)*inumlftn(irrep41)
     &                 + iposlftn(irrep1, irrep41) + (i1-1)*num4 + i4
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
      ELSEIF(TYPE.EQ.'1324')THEN
C
C            P Q R S -> P R Q S
C            1 2 3 4    1 3 2 4
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop1(irr2)*pop3(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop2(irr2)*pop4(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop2(irrep2r)*pop4(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop3(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep24 = dirprd(irrep2, irrep4)
                   irrep13 = dirprd(irrep1, irrep3)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
                 iold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + iposlfto(irrep2, irrep12)
                 inew = iabsoffn(irrep24)+iposrhtn(irrep4, irrep24)
     &                 + iposlftn(irrep3, irrep13)
                 nlfto = inumlfto(irrep12)
                 nlftn = inumlftn(irrep13)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 wout(inew
     &                 + ((i4-1)*num2+i2-1)*nlftn
     &                  + (i3-1)*num1 + i1)
     &          =win(iold
     &                 + ((i4-1)*num3+i3-1)*nlfto
     &                 + (i2-1)*num1 + i1)
C                ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
C    &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
C    &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
C                ioffnew = iabsoffn(irrep24)+iposrhtn(irrep4, irrep24)
C    &                 + ((i4-1)*num2+i2-1)*inumlftn(irrep13)
C    &                 + iposlftn(irrep3, irrep13) + (i3-1)*num1 + i1
C                wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'2134')THEN
C
C            P Q R S -> Q P R S
C            1 2 3 4    2 1 3 4
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop2(irr2)*pop1(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop3(irr2)*pop4(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop2(irrep2L)*pop1(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep34 = dirprd(irrep3, irrep4)
                   irrep21 = dirprd(irrep2, irrep1)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep34)+iposrhtn(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlftn(irrep21)
     &                 + iposlftn(irrep1, irrep21) + (i1-1)*num2 + i2
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
      ELSEIF(TYPE.EQ.'3124')THEN
C
C            P Q R S -> R P Q S
C            1 2 3 4    3 1 2 4
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop3(irr2)*pop1(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop2(irr2)*pop4(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop2(irrep2r)*pop4(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop3(irrep2L)*pop1(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep31 = dirprd(irrep3, irrep1)
                   irrep24 = dirprd(irrep2, irrep4)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep24)+iposrhtn(irrep4, irrep24)
     &                 + ((i4-1)*num2+i2-1)*inumlftn(irrep31)
     &                 + iposlftn(irrep1, irrep31) + (i1-1)*num3 + i3
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'3214')THEN
C
C            P Q R S -> R Q P S
C            1 2 3 4    3 2 1 4
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop3(irr2)*pop2(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop1(irr2)*pop4(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop1(irrep2r)*pop4(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop3(irrep2L)*pop2(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep32 = dirprd(irrep3, irrep2)
                   irrep14 = dirprd(irrep1, irrep4)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep14)+iposrhtn(irrep4, irrep14)
     &                 + ((i4-1)*num1+i1-1)*inumlftn(irrep32)
     &                 + iposlftn(irrep2, irrep32) + (i2-1)*num3 + i3
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'2413')THEN
C
C            P Q R S -> Q S P R 
C            1 2 3 4    2 4 1 3 
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop2(irr2)*pop4(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop1(irr2)*pop3(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop1(irrep2r)*pop3(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop2(irrep2L)*pop4(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep24 = dirprd(irrep2, irrep4)
                   irrep13 = dirprd(irrep1, irrep3)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep13)+iposrhtn(irrep3, irrep13)
     &                 + ((i3-1)*num1+i1-1)*inumlftn(irrep24)
     &                 + iposlftn(irrep4, irrep24) + (i4-1)*num2 + i2
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'1342')THEN
C
C            P Q R S -> P R S Q 
C            1 2 3 4    1 3 4 2 
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop1(irr2)*pop3(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop4(irr2)*pop2(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop4(irrep2r)*pop2(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop3(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep42 = dirprd(irrep4, irrep2)
                   irrep13 = dirprd(irrep1, irrep3)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = iabsoffo(irrep34) + iposrhto(irrep4, irrep34)
     &                 + ((i4-1)*num3+i3-1)*inumlfto(irrep12)
     &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
                 ioffnew = iabsoffn(irrep42)+iposrhtn(irrep2, irrep42)
     &                 + ((i2-1)*num4+i4-1)*inumlftn(irrep13)
     &                 + iposlftn(irrep3, irrep13) + (i3-1)*num1 + i1
C                 write(6,*)'sstgensowz',  ioffnew, ioffold
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
      ELSEIF(TYPE.EQ.'2341')THEN
C
C            P Q R S -> Q R S P
C            1 2 3 4    2 3 4 1
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop2(irr2)*pop3(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop4(irr2)*pop1(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop4(irrep2r)*pop1(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop2(irrep2L)*pop3(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
            enddo
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep23 = dirprd(irrep2, irrep3)
                   irrep41 = dirprd(irrep1, irrep4)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
                   ilfto = inumlfto(irrep12)
                   ilftn = inumlftn(irrep23)
                   ioffoo = iabsoffo(irrep34) + iposrhto(irrep4,irrep34)
     &                     + iposlfto(irrep2, irrep12)
                   ioffnn = iabsoffn(irrep41)+iposrhtn(irrep1, irrep41)
     &                     + iposlftn(irrep3, irrep23)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = ioffoo
     &                 + ((i4-1)*num3+i3-1)*ilfto
     &                 +  (i2-1)*num1 + i1
                 ioffnew = ioffnn
     &                 + ((i1-1)*num4+i4-1)*ilftn
     &                 + (i3-1)*num2 + i2
C    &                 + iposlfto(irrep2, irrep12) + (i2-1)*num1 + i1
C
C                ioffnew = iabsoffn(irrep14)+iposrhtn(irrep4, irrep14)
C    &                 + ((i4-1)*num1+i1-1)*inumlftn(irrep23)
C    &                 + iposlftn(irrep3, irrep23) + (i3-1)*num2 + i2
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      ELSEIF(TYPE.EQ.'3241')THEN
C
C            P Q R S -> R Q S P
C            1 2 3 4    3 2 4 1
C
C COMPUTE SIZE OF PS AND QR DISTRIBUTIONS FOR EACH IRREP
C
       call izero(inumlftn,8)
       call izero(inumrhtn,8)
       call izero(inumlfto,8)
       call izero(inumrhto,8)
       do irrep = 1, nirrep
          do irr1 = 1, nirrep
             irr2 = dirprd(irr1, irrep)
             inumlftn(irrep) = inumlftn(irrep)+pop3(irr2)*pop2(irr1)
             inumrhtn(irrep) = inumrhtn(irrep)+pop4(irr2)*pop1(irr1)
             inumlfto(irrep) = inumlfto(irrep)+pop1(irr2)*pop2(irr1)
             inumrhto(irrep) = inumrhto(irrep)+pop3(irr2)*pop4(irr1)
          enddo
C          write(6,*)'inum', inumlftn(irrep), inumrhtn(irrep),
C     &                      inumlfto(irrep), inumrhto(irrep)
       enddo
C
C NOW COMPUTE ABSOLUTE OFFSETS TO BEGINNING OF EACH DPD IRREP
C
       iabsoffn(1) = 0
       iabsoffo(1) = 0
C       write(6,*)'irrepx', irrepx
       do irrrht = 1, nirrep
          irrlft = dirprd(irrrht, irrepx)
          iabsoffn(irrrht+1) = iabsoffn(irrrht) +
     &                    inumlftn(irrlft)*inumrhtn(irrrht)
          iabsoffo(irrrht+1) = iabsoffo(irrrht) +
     &                    inumlfto(irrlft)*inumrhto(irrrht)
           iposrhto(1, irrrht) = 0
           iposlfto(1, irrlft) = 0
           iposrhtn(1, irrrht) = 0
           iposlftn(1, irrlft) = 0
C             write(6,*)'iposlftn', iposlftn(1,1)
           do irrep1 = 1, nirrep
              irrep2r = dirprd(irrep1, irrrht)
              irrep2l = dirprd(irrep1, irrlft)
              iposrhtn(irrep1+1, irrrht) = iposrhtn(irrep1, irrrht) +
     &                                     pop4(irrep2r)*pop1(irrep1)*
     &                                     inumlftn(irrlft)
              iposlftn(irrep1+1, irrlft) = iposlftn(irrep1, irrlft) +
     &                                     pop3(irrep2L)*pop2(irrep1)
              iposrhto(irrep1+1, irrrht) = iposrhto(irrep1, irrrht) +
     &                                     pop3(irrep2r)*pop4(irrep1)*
     &                                     inumlfto(irrlft)
              iposlfto(irrep1+1, irrlft) = iposlfto(irrep1, irrlft) +
     &                                     pop1(irrep2L)*pop2(irrep1)
C            write(6,*)'iposlftn', irrep1,
C            irrlft,iposlftn(irrep1,irrlft)
            enddo
C            write(6,*)'iabsoff', iabsoffn(irrrht), iabsoffo(irrrht) 
        enddo
       do irreprht = 1, nirrep
          irreplft = dirprd(irreprht, irrepx)
C          write(6,*)'test here'
          do irrep4 = 1, nirrep
             irrep3 = dirprd(irreprht, irrep4)
             irrep34 = irreprht
             irrep12 = irreplft
             num3 = pop3(irrep3)
             num4 = pop4(irrep4)
                do irrep2 =1, nirrep
                   irrep1 =dirprd(irreplft, irrep2)
                   irrep32 = dirprd(irrep3, irrep2)
                   irrep41 = dirprd(irrep1, irrep4)
                   num2 = pop2(irrep2)
                   num1 = pop1(irrep1)
                   ilfto = inumlfto(irrep12)
                   ioffoo = iabsoffo(irrep34) + iposrhto(irrep4,irrep34)
     &                     + iposlfto(irrep2, irrep12)
                   ilftn = inumlftn(irrep32)
                   ioffnn = iabsoffn(irrep41)+iposrhtn(irrep1, irrep41)
     &                     +iposlftn(irrep2, irrep32)
!$acc kernels
!$acc loop independent
                   do i4 = 1, num4
!$acc loop independent
                   do i3 = 1, num3
!$acc loop independent
                   do i2 = 1, num2
!$acc loop independent
                   do i1 = 1, num1
                 ioffold = ioffoo
     &                 + ((i4-1)*num3+i3-1)*ilfto
     &                 +  (i2-1)*num1 + i1
                 ioffnew = ioffnn
     &                 + ((i1-1)*num4+i4-1)*ilftn
     &                 +  (i2-1)*num3 + i3
                 wout(ioffnew) = win(ioffold)
                    enddo
                    enddo
                    enddo
                    enddo
!$acc end kernels
             enddo
          enddo
        enddo
C
      else 
        write(6,*)'error in sstgen', type 
      endif
C

C        write(6,*)'norm in sstgen', type, dnrm2(nsize, win, 1), dnrm2(nsize, wout, 1)
C
      RETURN
      END
C  
C  
