      SUBROUTINE DOSODIIS2_S(ITER, VECINOUT, LENGTH, MAXDIM, SCR, MAXCOR,
     &                    RESMDIIS)
C
      IMPLICIT NONE
      INTEGER, intent(in) :: ITER, LENGTH, MAXDIM
      INTEGER*8, intent(in) :: MAXCOR
C
      REAL*4, intent(inout) :: RESMDIIS
      REAL*4, intent(inout) :: SCR(MAXCOR), VECINOUT(LENGTH)
C
      INTEGER :: LENVEC, NDIMR, IDIMR
      COMMON/DIISDAT1S/LENVEC,NDIMR,IDIMR
C
      INTEGER :: I, J, K, N, I0, I11, I12, I2, IEND, MAXDIM1, INDXF
      REAL*4 :: FACT, RESMITER, X
C
      INTEGER, external :: IDAMAX
      REAL*4, external :: SDOT
      REAL*4, external :: SNRM2
C
      INDXF(I,J,N)=I+(J-1)*N
C
      MAXDIM1 = MAXDIM + 2
      LENVEC = LENGTH
C
      I0  = 1
      I11 = I0  + MAXDIM1*MAXDIM1
      I12 = I11 + MAXDIM1*MAXDIM1
      I2  = I12 + MAXDIM1
      IEND = I2 + LENVEC
      RESMDIIS = 1.0
C     
      IF (MAXCOR<IEND) THEN
         WRITE(6,*) 'Insufficient memory in DOSODIIS'
         STOP
      ENDIF
C
C      write(6,*)' test diis before clmoio'
      IF (ITER==1) THEN
         call clmoiodiis_s(scr, lenvec, 1, maxdim*2, 'DSTMP2_S', 4096)
C      write(6,*)' test diis  clmoio'
         IDIMR = 1 
         NDIMR = 1
         OPEN(211, FILE='TMPAMPTUDE_S', STATUS='UNKNOWN',
     &        FORM='UNFORMATTED')
         REWIND 211
         WRITE(211) (VECINOUT(I), I = 1, LENVEC)
         RETURN
      ENDIF
C      write(6,*)' test diis after clmoio'
C
      REWIND 211
      READ(211) (SCR(I), I = 1, LENVEC)
      call putdiis_s(scr, 1, 1, idimr+maxdim, 'DSTMP2_S')
C
      IF(NDIMR==1) WRITE(211) (VECINOUT(I), I = 1, LENVEC)
C
      CALL SAXPY(LENVEC,-1.0,SCR(1),1,VECINOUT,1)
      call putdiis_s(vecinout, 1, 1, idimr, 'DSTMP2_S')
      RESMITER = VECINOUT(IDAMAX(LENVEC, VECINOUT, 1))
      DO I = 1, NDIMR
         call getdiis_s(scr(i2),1, 1, i, 'DSTMP2_S')
         DO J = 1, I
         call getdiis_s(vecinout,1, 1, J, 'DSTMP2_S')
            X = SDOT(LENVEC, SCR(I2), 1, VECINOUT, 1)
            SCR(INDXF(I,J,MAXDIM))=X
            SCR(INDXF(J,I,MAXDIM))=X
         ENDDO
      ENDDO
C
      CALL DODIIS_S(SCR, SCR(I11), SCR(I12), MAXDIM, NDIMR)
C
      VECINOUT = 0.0
      DO I = 1, NDIMR
         FACT = SCR(I12+I-1)
         call getdiis_s(scr(I2),1, 1, I, 'DSTMP2_S')
         CALL SAXPY(LENVEC, FACT, SCR(I2), 1, VECINOUT, 1)
      ENDDO
      RESMDIIS = VECINOUT(IDAMAX(LENVEC, VECINOUT, 1))
C
      WRITE(6,1001) ITER
      WRITE(6,1002) RESMITER
      WRITE(6,1003) RESMDIIS
C
      DO I = 1, NDIMR
         FACT = SCR(I12+I-1)
         call getdiis_s(scr(I2),1, 1, I+maxdim, 'DSTMP2_S')
         CALL SAXPY(LENVEC, FACT, SCR(I2), 1, VECINOUT, 1)
      ENDDO
C
      REWIND 211
      WRITE(211) (VECINOUT(I), I = 1, LENVEC)
C
      NDIMR=MIN(MAXDIM,NDIMR+1)
      IDIMR=IDIMR+1
      IF(IDIMR.GT.MAXDIM) IDIMR=1
C      write(6,*)'data in doso', LENVEC, NDIMR, IDIMR
C
1001  FORMAT(T3,'Convergence information after ',I5,' iterations: ')
1002  FORMAT(T3,'Largest element of residual vector : ',E15.8,'.')
1003  FORMAT(T3,'Largest element of DIIS residual   : ',E15.8,'.')
C
      RETURN
      END      

