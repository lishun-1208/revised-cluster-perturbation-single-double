      SUBROUTINE XGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
      CALL XGEMMD(TRANSA,TRANSB,M,N,K,ALPHA,A,MAX(1,LDA),B,MAX(1,LDB),
     $                   BETA, C, MAX(1,LDC) )
      RETURN
      END
