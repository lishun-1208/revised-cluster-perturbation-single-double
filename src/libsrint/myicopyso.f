
      SUBROUTINE MYICOPYSO(I1,I2,LEN)
C
C THIS ROUTINE COPIES THE FIRST LEN ELEMENTS OF INTEGER
C   VECTOR I1 INTO VECTOR I2.
      use mod_iop
      DOUBLE PRECISION, intent(in),managed:: I1(*)
      DOUBLE PRECISION, intent(out),managed:: I2(*)
      INTEGER, intent(in) ::LEN
C
      IF(PREC==8.and.igpu==0) THEN 
         CALL DCOPY(LEN, I1, 1, I2, 1)
      ELSEIF(PREC==4.and.igpu==0) THEN 
         CALL SCOPY(LEN, I1, 1, I2, 1) 
      ELSEIF(PREC==8.and.igpu==1) THEN 
         CALL DCOPYG(LEN, I1, 1, I2, 1)
      ELSEIF(PREC==4.and.igpu==1) THEN 
         CALL SCOPYG(LEN, I1, 1, I2, 1) 
      ENDIF
      RETURN
      END
!
      subroutine dcopyg(nsize, iin, ioffin, iout, ioffout)
      use cublas, only:cublasdcopy, cublasscopy 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*8, intent(in),managed :: iin(*), iout(*) 
      integer*4 ioffin4, ioffout4, nsize4 
C
      ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
      call cublasdcopy(nsize4, iin, ioffin4, iout, ioffout4)  
      return
      end
       
      subroutine scopyg(nsize, iin, ioffin, iout, ioffout)
      use cublas, only:cublasdcopy, cublasscopy 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*8, intent(in),managed :: iin(*), iout(*) 
      integer*4 ioffin4, ioffout4, nsize4 
C
      ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
      call cublasscopy(nsize4, iin, ioffin4, iout, ioffout4)  
      return
      end
       
