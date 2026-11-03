      SUBROUTINE XGEMMSO(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $                   BETA, C, LDC )
      use mod_iop
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION, managed :: A( LDA, * ), B( LDB, * ), C( LDC, * )
      REAL*4 ALPHA4, BETA4 
      if(prec==8.and.igpu==0) then 
      CALL XGEMMD(TRANSA,TRANSB,M,N,K,ALPHA,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA, C, MAX(1,LDC) )
      elseif(prec==4.and.igpu==0) then 
      ALPHA4 = ALPHA; BETA4 = BETA
      CALL XGEMMS(TRANSA,TRANSB,M,N,K,ALPHA4,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA4, C, MAX(1,LDC) )
      elseif(prec==8.and.igpu==1) then 
      CALL XGEMMDG(TRANSA,TRANSB,M,N,K,ALPHA,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA, C, MAX(1,LDC) )
      elseif(prec==4.and.igpu==1) then 
      CALL XGEMMSG(TRANSA,TRANSB,M,N,K,ALPHA,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA, C, MAX(1,LDC) )
      endif
      RETURN
      END
      SUBROUTINE XGEMMS(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL*4   ALPHA, BETA
*     .. Array Arguments ..
      REAL*4   A( LDA, * ), B( LDB, * ), C( LDC, * )

      CALL SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA, C, MAX(1,LDC) )
      
      RETURN
      END
      SUBROUTINE XGEMMD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )

      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA, C, MAX(1,LDC) )
      
      RETURN
      END
      SUBROUTINE XGEMMSG(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $                   BETA, C, LDC )
      use cublas 
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER*8          M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      real*4,device :: A( LDA, * ), B( LDB, * ), C( LDC, * )
C
      real*4 alpha4, beta4
      INTEGER*4 :: M4, N4, K4, LDA4, LDB4, LDC4, ONE 
      M4 = M; N4 = N; K4 = K 
      ONE = 1 
      LDA4 = MAX(LDA, 1); LDB4 = MAX(LDB, 1); LDC4 = MAX(LDC,1)
      alpha4 = alpha; beta4 = beta
C
      CALL CUBLASSGEMM(TRANSA,TRANSB,M4,N4,K4,ALPHA4,A,LDA4,
     &                   B,LDB4,BETA4, C, LDC4 )
      
      RETURN
      END
      SUBROUTINE XGEMMDG(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $                   BETA, C, LDC )
      use cublas 
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER*8          M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION,device :: A( LDA, * ), B( LDB, * ), C( LDC, * )
C
      INTEGER*4 :: M4, N4, K4, LDA4, LDB4, LDC4, ONE 
      M4 = M; N4 = N; K4 = K 
      ONE = 1 
      LDA4 = MAX(LDA, 1); LDB4 = MAX(LDB, 1); LDC4 = MAX(LDC,1)
C
      CALL CUBLASDGEMM(TRANSA,TRANSB,M4,N4,K4,ALPHA,A,LDA4,
     &                   B,LDB4,BETA, C, LDC4 )
      
      RETURN
      END
