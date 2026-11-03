      SUBROUTINE DOSODIIS2_D(ITER, VECINOUT, LENGTH, MAXDIM, SCR, MAXCOR,
     &                    RESMDIIS)
C
      IMPLICIT NONE
      INTEGER, intent(in) :: ITER, LENGTH, MAXDIM
      INTEGER*8, intent(in) :: MAXCOR
C
      DOUBLE PRECISION, intent(inout) :: RESMDIIS
      DOUBLE PRECISION, intent(inout) :: SCR(MAXCOR), VECINOUT(LENGTH)
C
      INTEGER :: LENVEC, NDIMR, IDIMR
      COMMON/DIISDAT1/LENVEC,NDIMR,IDIMR
C
      INTEGER :: I, J, K, N, I0, I11, I12, I2, IEND, MAXDIM1, INDXF
      DOUBLE PRECISION :: FACT, RESMITER, X
C
      INTEGER, external :: IDAMAX
      DOUBLE PRECISION, external :: DDOT
      DOUBLE PRECISION, external :: DNRM2
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
      RESMDIIS = 1.D0
C      write(6,*)'I in dosodiis', i11, i12, i2, iend
C     
      IF (MAXCOR<IEND) THEN
         WRITE(6,*) 'Insufficient memory in DOSODIIS'
         STOP
      ENDIF
C
C      write(6,*)'begin dosodiis'
      IF (ITER==1) THEN
         call clmoiodiis(scr, lenvec, 1, maxdim*2, 'DIISTMP2', 4096)
C         OPEN(20, FILE='DIIS.tmp', STATUS='UNKNOWN', ACCESS='DIRECT',
C     &        FORM='UNFORMATTED', RECL=8*LENGTH)
         IDIMR = 1 
         NDIMR = 1
         OPEN(21, FILE='TMPAMPTUDE', STATUS='UNKNOWN',
     &        FORM='UNFORMATTED')
         REWIND 21
         WRITE(21) (VECINOUT(I), I = 1, LENVEC)
         RETURN
      ENDIF
C
      REWIND 21
      READ(21) (SCR(I), I = 1, LENVEC)
C      write(6,*)'begin dosodiis'
C      write(6,*)'rec and ndimr',idimr+maxdim, ndimr
      call putdiis(scr, 1, 1, idimr+maxdim, 'DIISTMP2')
C      WRITE(20, REC=IDIMR+MAXDIM) (SCR(I), I = 1, LENVEC)     
C
      IF(NDIMR==1) WRITE(21) (VECINOUT(I), I = 1, LENVEC)
C
C      write(6,*)'begin dosodiis'
      CALL DAXPY(LENVEC,-1.D0,SCR(1),1,VECINOUT,1)
C      WRITE(20, REC=IDIMR) (VECINOUT(I), I = 1, LENVEC)     
      call putdiis(vecinout, 1, 1, idimr, 'DIISTMP2')
      RESMITER = VECINOUT(IDAMAX(LENVEC, VECINOUT, 1))
C      write(6,*)'begin dosodiis'
      DO I = 1, NDIMR
         call getdiis(scr(i2),1, 1, i, 'DIISTMP2')
C         READ(20, REC=I) (SCR(K), K = I2, IEND-1)
         DO J = 1, I
         call getdiis(vecinout,1, 1, J, 'DIISTMP2')
C            READ(20, REC=J) (VECINOUT(K), K = 1, LENVEC)
            X = DDOT(LENVEC, SCR(I2), 1, VECINOUT, 1)
C            write(6,'(A, I3,I3," ", E15.8)')'x in dosodiis', i, j, x
            SCR(INDXF(I,J,MAXDIM))=X
            SCR(INDXF(J,I,MAXDIM))=X
         ENDDO
      ENDDO
C
      CALL DODIIS(SCR, SCR(I11), SCR(I12), MAXDIM, NDIMR)
C
      VECINOUT = 0.D0
      DO I = 1, NDIMR
         FACT = SCR(I12+I-1)
C        WRITE(6,*) 'I COEF', I, FACT
C         READ(20, REC=I) (SCR(K), K = I2, IEND-1)
         call getdiis(scr(I2),1, 1, I, 'DIISTMP2')
         CALL DAXPY(LENVEC, FACT, SCR(I2), 1, VECINOUT, 1)
      ENDDO
      RESMDIIS = VECINOUT(IDAMAX(LENVEC, VECINOUT, 1))
C
      WRITE(6,1001) ITER
      WRITE(6,1002) RESMITER
      WRITE(6,1003) RESMDIIS
C
      DO I = 1, NDIMR
         FACT = SCR(I12+I-1)
         call getdiis(scr(I2),1, 1, I+maxdim, 'DIISTMP2')
C         READ(20, REC=I+MAXDIM) (SCR(K), K = I2, IEND-1)
         CALL DAXPY(LENVEC, FACT, SCR(I2), 1, VECINOUT, 1)
      ENDDO
C
      REWIND 21
      WRITE(21) (VECINOUT(I), I = 1, LENVEC)
C
      NDIMR=MIN(MAXDIM,NDIMR+1)
      IDIMR=IDIMR+1
      IF(IDIMR.GT.MAXDIM) IDIMR=1
C      write(6,*)'data in doso', LENVEC, NDIMR, IDIMR
C
C      write(6,*)'dosodiis', idimr, maxdim, ndimr
1001  FORMAT(T3,'Convergence information after ',I5,' iterations: ')
1002  FORMAT(T3,'Largest element of residual vector : ',E15.8,'.')
1003  FORMAT(T3,'Largest element of DIIS residual   : ',E15.8,'.')
C
      RETURN
      END      

