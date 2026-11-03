      SUBROUTINE TRANSPSO(A,B,NUM,DIS)
C
C    THIS ROUTINE PERFORMS A GENERAL MATRIX  TRANSPOSITION 
C
C      A(Q,P) --> B(P,Q) 
C
C    INPUT:
C             A  ......   MATRIX TO BE TRANSPOSED
C             NUM .....   LENGTH OF SECOND INDEX OF A
C                         (NUMBER OF DISTRIBUTIONS IN A)
C             DIS .....   LENGTH OF FIRST INDEX OF A
C                         (DISTRIBUTION LENGTH IN A)
C
C    OUTPUT:
C             B .......   TRANSPOSED MATRIX
C
CEND
C
C  CODED JULY/90  JG
C
      use omp_lib
      use mod_iop
      IMPLICIT NONE 
      INTEGER DIS, NUM 
      real*4, intent(in),managed ::  A(*)
      real*4, intent(out),managed :: B(*)
C     this shows much better performance MEH August 2006
      if(prec==8.and.igpu==0) then 
        call TRANSPSOd(A,B,NUM,DIS)
      elseif(prec==4.and.igpu==0) then 
        call TRANSPSOs(A,B,NUM,DIS)
      elseif(prec==8.and.igpu==1) then 
        call TRANSPSOdg(A,B,NUM,DIS)
      elseif(prec==4.and.igpu==1) then 
        call TRANSPSOsg(A,B,NUM,DIS)
      endif
      RETURN
      END
      SUBROUTINE TRANSPSOS(A,B,NUM,DIS)
C
C    THIS ROUTINE PERFORMS A GENERAL MATRIX  TRANSPOSITION 
C
C      A(Q,P) --> B(P,Q) 
C
C    INPUT:
C             A  ......   MATRIX TO BE TRANSPOSED
C             NUM .....   LENGTH OF SECOND INDEX OF A
C                         (NUMBER OF DISTRIBUTIONS IN A)
C             DIS .....   LENGTH OF FIRST INDEX OF A
C                         (DISTRIBUTION LENGTH IN A)
C
C    OUTPUT:
C             B .......   TRANSPOSED MATRIX
C
CEND
C
C  CODED JULY/90  JG
C
      use omp_lib
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS 
      real*4, intent(in) ::  A(*)
      real*4, intent(out) :: B(*)
C          A(DIS,NUM),B(NUM,DIS)
C
COLD      DO 10 IJ=1,NUM
COLD      DO 10 IAB=1,DIS
COLD       B(IJ,IAB)=A(IAB,IJ)
COLD10    CONTINUE
C

C     this shows much better performance MEH August 2006
!$omp parallel do 
      DO IAB=1,DIS
      DO IJ=1,NUM
C      DO 10 IAB=1,DIS
       B((IAB-1)*NUM+IJ)=A((IJ-1)*DIS+IAB)
C       B(IJ,IAB)=A(IAB,IJ)
      ENDDO
      ENDDO
!$omp end parallel do 
      RETURN
      END
      SUBROUTINE TRANSPSOD(A,B,NUM,DIS)
C
C    THIS ROUTINE PERFORMS A GENERAL MATRIX  TRANSPOSITION 
C
C      A(Q,P) --> B(P,Q) 
C
C    INPUT:
C             A  ......   MATRIX TO BE TRANSPOSED
C             NUM .....   LENGTH OF SECOND INDEX OF A
C                         (NUMBER OF DISTRIBUTIONS IN A)
C             DIS .....   LENGTH OF FIRST INDEX OF A
C                         (DISTRIBUTION LENGTH IN A)
C
C    OUTPUT:
C             B .......   TRANSPOSED MATRIX
C
CEND
C
C  CODED JULY/90  JG
C
      use omp_lib
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS 
      DOUBLE PRECISION, intent(in) ::  A(*)
      DOUBLE PRECISION, intent(out) :: B(*)
C          A(DIS,NUM),B(NUM,DIS)
C
COLD      DO 10 IJ=1,NUM
COLD      DO 10 IAB=1,DIS
COLD       B(IJ,IAB)=A(IAB,IJ)
COLD10    CONTINUE
C

C     this shows much better performance MEH August 2006
!$omp parallel do 
      DO IAB=1,DIS
      DO IJ=1,NUM
C      DO 10 IAB=1,DIS
       B((IAB-1)*NUM+IJ)=A((IJ-1)*DIS+IAB)
C       B(IJ,IAB)=A(IAB,IJ)
      ENDDO
      ENDDO
!$omp end parallel do 
      RETURN
      END
      SUBROUTINE TRANSPSOSG(A,B,NUM,DIS)
C
C    THIS ROUTINE PERFORMS A GENERAL MATRIX  TRANSPOSITION 
C
C      A(Q,P) --> B(P,Q) 
C
C    INPUT:
C             A  ......   MATRIX TO BE TRANSPOSED
C             NUM .....   LENGTH OF SECOND INDEX OF A
C                         (NUMBER OF DISTRIBUTIONS IN A)
C             DIS .....   LENGTH OF FIRST INDEX OF A
C                         (DISTRIBUTION LENGTH IN A)
C
C    OUTPUT:
C             B .......   TRANSPOSED MATRIX
C
CEND
C
C  CODED JULY/90  JG
C
      use omp_lib
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS 
      real*4, intent(in), device:: A(DIS,NUM)
      real*4, intent(out), device:: B(NUM,DIS)
C
COLD      DO 10 IJ=1,NUM
COLD      DO 10 IAB=1,DIS
COLD       B(IJ,IAB)=A(IAB,IJ)
COLD10    CONTINUE
C

C     this shows much better performance MEH August 2006
!$acc parallel loop
      DO IAB=1,DIS
      DO IJ=1,NUM
C      DO 10 IAB=1,DIS
       B(IJ,IAB)=A(IAB,IJ)
      ENDDO
      ENDDO
!$acc end parallel loop 
      RETURN
      END
      SUBROUTINE TRANSPSODG(A,B,NUM,DIS)
C
C    THIS ROUTINE PERFORMS A GENERAL MATRIX  TRANSPOSITION 
C
C      A(Q,P) --> B(P,Q) 
C
C    INPUT:
C             A  ......   MATRIX TO BE TRANSPOSED
C             NUM .....   LENGTH OF SECOND INDEX OF A
C                         (NUMBER OF DISTRIBUTIONS IN A)
C             DIS .....   LENGTH OF FIRST INDEX OF A
C                         (DISTRIBUTION LENGTH IN A)
C
C    OUTPUT:
C             B .......   TRANSPOSED MATRIX
C
CEND
C
C  CODED JULY/90  JG
C
      use omp_lib
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS 
      DOUBLE PRECISION, intent(in), device:: A(DIS,NUM)
      DOUBLE PRECISION, intent(out), device:: B(NUM,DIS)
C
COLD      DO 10 IJ=1,NUM
COLD      DO 10 IAB=1,DIS
COLD       B(IJ,IAB)=A(IAB,IJ)
COLD10    CONTINUE
C

C     this shows much better performance MEH August 2006
!$acc parallel loop
      DO IAB=1,DIS
      DO IJ=1,NUM
C      DO 10 IAB=1,DIS
       B(IJ,IAB)=A(IAB,IJ)
      ENDDO
      ENDDO
!$acc end parallel loop 
      RETURN
      END
