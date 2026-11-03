C
      SUBROUTINE MINVR(AB,N,LDA,SCRATCH,DET,EPS,M,MODE)
C
C     A subroutine that calculates the determinant and inverse of
C          a matrix, as well as solving systems of linear equations.
C     Martin J. McBride.  11/25/85.
C     General Electric CRD, Information System Operation.
C
C  THE COMMENTS ABOVE REFER TO AN OLD "SCIPORT" VERSION OF THIS ROUTINE,
C  PRESUMABLY WRITTEN BY MR. MC BRIDE.  UNFORTUNATELY, IT SEEMS THAT 
C  HE IS A BAD PROGRAMMER AND THIS ROUTINE DID NOT WORK.  MINV IS NOW
C  A SIMPLE INTERFACE TO LINPACK ROUTINES WHICH *WORK*
C
CEND
      IMPLICIT NONE
      integer, intent(in) :: N, LDA,  M, MODE
      real*4, intent(in) :: EPS, DET(2)
      real*4, intent(inout) :: AB(LDA, LDA+1), SCRATCH(*)
C
      integer ioff
C
      IOFF=1+N
C      write(6,*)'test minvr'
C      write(6,*)'N', N, LDA, M, MODE
      CALL SGECO(AB,LDA,N,SCRATCH,DET,SCRATCH(IOFF))
C      write(6,*)'test minvr'
      IF(M.EQ.1)CALL SGESL(AB,LDA,N,SCRATCH,AB(1,N+1),0)
      IF(M.GT.1)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@MINV-F, Only one system can be solved at a time.')
C       CALL ERREX
      ENDIF
      IF(M.EQ.0)CALL SGEDI(AB,LDA,N,SCRATCH,DET,SCRATCH(IOFF),11)
C
      RETURN
      END
