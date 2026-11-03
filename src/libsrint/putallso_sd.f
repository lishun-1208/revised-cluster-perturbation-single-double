      SUBROUTINE PUTALLSO_SD(Z,list,filename)
C
C THIS ROUTINE WRITES OUT ALL IRREPS OF A SYMMETRY PACKED DISTRIBUTION.
C
C      Z  - VECTOR HOLDING THE LIST.
C      filename - which file to put vector Z
C
CEND
      use mod_iop
C
      integer, intent(in) :: list
      DOUBLE PRECISION, intent(inout):: Z(*)
      character(len=8), intent(in) :: filename
C
      IF(PREC==8) THEN
         CALL PUTALLSO_D(Z, LIST, FILENAME)
      ELSEIF(PREC==4) THEN
         CALL PUTALLSO_S(Z, LIST, FILENAME)
      ENDIF
C
      RETURN
      END
