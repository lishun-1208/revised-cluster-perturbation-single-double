      SUBROUTINE GETLSTSO(ZLIST,NFIRST,NBLOCK,
     &                     ISPIN,ILIST,filename)
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
C     NFIRST- THE NUMBER OF THE FIRST DISTRIBUTION WHICH WILL BE
C             WRITTEN TO OR READ FROM DISK.
C     NBLOCK- THE NUMBER OF DISTRIBUTIONS WHICH WILL BE HANDLED.
C             T AMPLITUDE I/O AND 2 FOR INTEGRAL I/O.
C     ISPIN - THE SYMMETRY OF THE RIGHT SIDE (IRREPRS) OF WPQRS 
C     ILIST - THE LIST of the saved Wpqrs with defined spincase. 
C
CEND
      use mod_symm
      use mod_orbit
      use mod_size
      IMPLICIT NONE 
      integer, intent(in)::nfirst, nblock, ispin, ilist 
      character(len=8), intent(in) :: filename
     
      double precision, intent(inout) :: ZLIST(*)
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      DOUBLE PRECISION, ALLOCATABLE :: TMP(:)
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
C
      integer(kind=8) :: nwords, nstart, nend, nwrdf, nwrdf2
      integer*8 test1, test2
      integer irrepx, i, j, ifile, listrec
      integer recstart,  recend, istat
      integer iout, irecord, istart
      character*200 filepath

      call get_f_path(filepath)

C
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      OPEN(UNIT=413, FILE=trim(filepath)//trim(FILENAME2),
     &           FORM='UNFORMATTED')
      rewind 413
      read(413) IRREPX
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
C      do j = 1, 20
C      do i = 1, 10 
C         write(6,*)'rectmp(i,j)', i, j, rectmp(i,j)
C      enddo
C      enddo
C      write(6,*)'end read ', filename2, ' in getlst', iprcln
      OPEN(UNIT=412, FILE=trim(filepath)//trim(FILENAME), 
     &            FORM='UNFORMATTED', access='direct',
     &            recl=iprcln)       
C
C DEAL WITH RIDICULOUS LIST NUMBERS
C
      ifile = 412 
      listrec = rectmp(ispin, ilist)
C define the started rec and ended rec
        NWORDS=NBLOCK*MOIONUML(ISPIN,ILIST)
      IF(NWORDS<=0) THEN
         close (412)
         RETURN
      ENDIF
        nstart = (nfirst-1)*moionuml(ispin, ilist)
        nend = (nfirst+nblock-1)*moionuml(ispin, ilist)
          RECSTART=(NSTART-1)/IPRCWD
          NWRDF = NSTART - iprcwd*recstart
          RECSTART = RECSTART + LISTREC + 1
C
          RECEND = (NEND-1)/IPRCWD
          nwrdf2 = nend - iprcwd*recend
          RECEND = RECEND + LISTREC + 1
C
          ALLOCATE(TMP(ICHCSZ), STAT=ISTAT)
C read the first line
C          write(6,100)'moionum',moionuml(ispin, ilist), nfirst-1,
C     &          nstart, icrsiz
C          write(6,100)'moionum',moionuml(ispin, ilist)
C          write(6,100)'nwords', nwords, nfirst, nstart, nend, recstart
C          write(6,100)'nwords', recend, nend, nwrdf2
100       FORMAT(A, I12, I12, I12, I12, I12)
          IF(RECEND>RECSTART) THEN 
C               write(6,*)'irecord', irecord
            CALL RDDIRSO(IFILE, RECSTART, TMP, ICHCSZ)   
            CALL MYICOPY(TMP(NWRDF+1),ZLIST, ICHCSZ-NWRDF)
C read from the second to the penul-timate line 
            DO IOUT = 1, RECEND-RECSTART-1
               IRECORD = RECSTART + IOUT 
C               write(6,*)'irecord', irecord
               ISTART = (iout-1)*ICHCSZ + ICHCSZ - NWRDF + 1
               CALL RDDIRSO(IFILE,IRECORD,ZLIST(istart),ICHCSZ)
            ENDDO
C read the last line 
            CALL RDDIRSO(IFILE, RECEND, TMP, ICHCSZ)   
            istart = nwords - nwrdf2 + 1
C            write(6,*)'istart', istart
            CALL MYICOPY(TMP, ZLIST(istart), NWRDF2)
          ELSEIF(RECEND==RECSTART) THEN 
C               write(6,*)'irecord',recstart, nwrdf 
            CALL RDDIRSO(IFILE, RECSTART, TMP, ICHCSZ) 
            CALL MYICOPY(TMP(NWRDF+1), ZLIST, NWRDF2-NWRDF) 
          ENDIF
          deallocate(TMP)
      close(412)
C write(6,*)'========end of getlst======='
      RETURN
      END
