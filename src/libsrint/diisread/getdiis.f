      SUBROUTINE GETDIIS(ZLIST,NFIRST,NBLOCK,ILIST,filename)
CEND
      IMPLICIT NONE
      integer, intent(in)::nfirst, nblock, ilist
      character(len=8), intent(in) :: filename
C
      double precision, intent(inout) :: ZLIST(*) 
C     
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(1000), MOIONUMR(1000), RECTMP(1000)
C
      double precision, allocatable :: TMP(:)
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
C
      integer*8 nwords, nstart, nend, nwrdf, nwrdf2
      integer*8 numl
      integer numr, dsize, nlist 
      integer irrepx, i, j, ifile, listrec
      integer recstart, recend, istat
      integer iout, irecord, istart
      integer*8 tmpi1, tmpi2
C
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      filename2=trim(filename)//tmpc
C      write(6,*)'begin read',' ', filename2, ' in putlsttmp'
      OPEN(UNIT=413, FILE=FILENAME2, FORM='UNFORMATTED',status='old')
      rewind 413
      read(413) nlist 
      read(413) dsize
      read(413) numl
      read(413) numr
      read(413) (rectmp(i),i=1,nlist)
      close(413)
C
C DEAL WITH RIDICULOUS LIST NUMBERS
C
      OPEN(UNIT=415, FILE=FILENAME, FORM='UNFORMATTED', access='direct',
     &            recl=dsize*8)
      ifile = 415 
         listrec = rectmp(ilist)
C define the started rec and ended rec , one should note that recstart
C and recend are both integer
        tmpi1 = nblock; tmpi2 = numl 
        nwords = tmpi1*tmpi2
C
        tmpi1 = nfirst-1 
        nstart = tmpi1*tmpi2
C
        tmpi1 = nfirst+nblock-1
        nend = tmpi1*tmpi2
C
        recstart=(nstart-1)/dsize
        tmpi1 = dsize; tmpi2 = recstart 
        nwrdf = nstart - tmpi1*tmpi2
        recstart = recstart + listrec + 1
C 
        recend = (nend-1)/dsize
        tmpi2 = recend
        nwrdf2 = nend - tmpi1*tmpi2
        recend = recend + listrec + 1
C
        allocate(tmp(dsize), stat = istat) 
C write the first line
        if(recend>recstart) then 
           call rddirso(ifile, recstart, tmp, dsize)
           call myicopy(TMP(NWRDF+1), ZLIST, dsize-NWRDF)
C write the second to the penul-timate line
           do iout = 1, recend-recstart-1 
              irecord = recstart + iout
              istart = (iout-1)*dsize + dsize - NWRDF + 1 
              CALL rddirso(IFILE,IRECORD,ZLIST(istart),dsize)
           enddo 
C write the last line 
        call rddirso(ifile, recend, tmp, dsize)
        istart = nwords - nwrdf2 + 1
        call myicopy(tmp, Zlist(istart), nwrdf2)
C
        elseif(recend==recstart) then             
          call rddirso(ifile, recstart, tmp, dsize)
          call myicopy(tmp(nwrdf+1), zlist, nwrdf2-nwrdf)
        endif
        deallocate(TMP)
C       ENDIF
C      write(6,*)'========end of putlst======='
      close(415)
300   RETURN
      END
