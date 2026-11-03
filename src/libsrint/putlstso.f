      SUBROUTINE PUTLSTSO(ZLIST,NFIRST,NBLOCK,
     &                           ISPIN,ILIST,filename)
C
C THIS SUBROUTINE RETRIEVES NBLOCK DISTRIBUTIONS OF A PARTICULAR INTEGRA
C  LIST, BEGINNING WITH DISTRIBUTION NFIRST.  THE VALUES ARE RETURNED IN
C  VECTOR ZLIST.
C
C   OUTPUT:
C
C         ZLIST - THE TARGET VECTOR USED IN I/O OPERATIONS.
C
C   INPUT PARAMETERS:
C
C         NFIRST- THE NUMBER OF THE FIRST DISTRIBUTION WHICH WILL BE
C                 WRITTEN TO OR READ FROM DISK.
C         NBLOCK- THE NUMBER OF DISTRIBUTIONS WHICH WILL BE HANDLED.
C         USECIN- THE CACHE BUFFER WHICH IS TO BE USED FOR I/O OPERATION
C                 BEST PERFORMANCE IS ACHIEVED IF THIS IS SET TO 1 FOR
C                 T AMPLITUDE I/O AND 2 FOR INTEGRAL I/O.
C         ISPIN - THE LEFT-HAND INDEX OF THE LIST [MOIO(ISPIN,ITYPE)]. ispin=symQ, for I(P,Q)
C         ILIST - THE RIGHT-HAND INDEX OF THE LIST [MOIO(ISPIN,ITYPE)]. ilist=1
C
CEND
      use mod_symm
      use mod_size
      IMPLICIT NONE
      integer, intent(in)::nfirst, nblock, ispin, ilist 
      character(len=8), intent(in) :: filename
C
      double precision, intent(inout) :: ZLIST(*) 
C     
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      double precision, allocatable :: TMP(:)
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
C
      integer*8 nwords, nstart, nend, nwrdf, nwrdf2
      integer irrepx, i, j, ifile, listrec
      integer recstart, recend, istat
      integer iout, irecord, istart
      integer*8 tmpi1, tmpi2
      character*200 filepath

      call get_f_path(filepath)
C
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      filename2=trim(filename)//tmpc
C      write(6,*)'begin read',' ', filename2, ' in putlsttmp'
      OPEN(UNIT=413, FILE=trim(filepath)//trim(FILENAME2), 
     &          FORM='UNFORMATTED',status='old')
      rewind 413
      read(413) IRREPX
C      write(6,*)'test'
      read(413) (popp(i), i=1,8)
      read(413) (popq(i), i=1,8)
      read(413) (popr(i), i=1,8)
      read(413) (pops(i), i=1,8)
      read(413) (ioppq(i), i = 1,8)
      read(413) (ioprs(i), i = 1,8)
      read(413) (irrep0(i), i = 1,8)
      read(413) (irrep0(i), i = 9,16)
      read(413) ((moionuml(i,j),i=1,10),j=1,20)
      read(413) ((moionumr(i,j),i=1,10),j=1,20)
      read(413) ((rectmp(i,j),i=1,10),j=1,20)
      close(413)
C      write(6,*)'end read',' ', trim(filename2), ' in putlsttmp'
C      goto 300
C
C DEAL WITH RIDICULOUS LIST NUMBERS
C
      OPEN(UNIT=415, FILE=trim(filepath)//trim(FILENAME), 
     &            FORM='UNFORMATTED', access='direct',
     &            recl=iprcln)
      ifile = 415 
         listrec = rectmp(ispin, ilist)
C define the started rec and ended rec , one should note that recstart
C and recend are both integer
        nwords = nblock*moionuml(ispin, ilist)
      IF(nwords<=0) THEN
         close(415)
         RETURN
      ENDIF
        nstart = (nfirst-1)*moionuml(ispin, ilist) 
        nend = (nfirst+nblock-1)*moionuml(ispin, ilist) 
        recstart=(nstart-1)/iprcwd
        nwrdf = nstart - iprcwd*recstart
        recstart = recstart + listrec + 1
C 
        recend = (nend-1)/iprcwd
        nwrdf2 = nend - iprcwd*recend
        recend = recend + listrec + 1
C
        allocate(tmp(ichcsz), stat = istat) 
C write the first line
        if(recend>recstart) then 
           call rddirso(ifile, recstart, tmp, ichcsz)
           call myicopy(ZLIST, TMP(NWRDF+1), ICHCSZ-NWRDF)
           call wrdirso(ifile, recstart, tmp, ichcsz)
C write the second to the penul-timate line
           do iout = 1, recend-recstart-1 
              irecord = recstart + iout
              istart = (iout-1)*ICHCSZ + ICHCSZ - NWRDF + 1 
              CALL wrdirso(IFILE,IRECORD,ZLIST(istart),ICHCSZ)
           enddo 
C write the last line 
        call rddirso(ifile, recend, tmp, ichcsz)
        istart = nwords - nwrdf2 + 1
        call myicopy(Zlist(istart), tmp, nwrdf2)
        CALL wrdirso(IFILE,recend,tmp,ICHCSZ)
C
        elseif(recend==recstart) then             
          call rddirso(ifile, recstart, tmp, ichcsz)
          call myicopy(zlist, tmp(nwrdf+1), nwrdf2-nwrdf)
          CALL wrdirso(IFILE,recstart,tmp,ICHCSZ)
        endif
        deallocate(TMP)
C       ENDIF
C      write(6,*)'========end of putlst======='
      close(415)
300   RETURN
      END
