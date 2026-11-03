      SUBROUTINE PUTDIIS(ZLIST,NFIRST,NBLOCK,ILIST,filename)
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
      integer numr, dsize, iopsize
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
!      write(6,*)'begin read',' ', filename2, ' in putlsttmp'
      OPEN(UNIT=413, FILE=trim(filepath)//trim(FILENAME2), 
     &                FORM='UNFORMATTED',status='old')
      rewind 413
      read(413) iopsize 
      read(413) dsize
      read(413) numl
      read(413) numr
      read(413) (rectmp(i),i=1,iopsize)
      close(413)
C
C DEAL WITH RIDICULOUS LIST NUMBERS
C
      OPEN(UNIT=415, FILE=trim(filepath)//trim(FILENAME), 
     &            FORM='UNFORMATTED', access='direct',
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
           call myicopy(ZLIST, TMP(NWRDF+1), dsize-NWRDF)
           call wrdirso(ifile, recstart, tmp, dsize)
C write the second to the penul-timate line
           do iout = 1, recend-recstart-1 
              irecord = recstart + iout
              istart = (iout-1)*dsize + dsize - NWRDF + 1 
              CALL wrdirso(IFILE,IRECORD,ZLIST(istart),dsize)
           enddo 
C write the last line 
C           write(6,*)'rddir1', iout, recend
        call rddirso(ifile, recend, tmp, dsize)
C           write(6,*)'rddir1', iout, recend
        istart = nwords - nwrdf2 + 1
        call myicopy(Zlist(istart), tmp, nwrdf2)
C           write(6,*)'rddir1', iout, recend
        CALL WRDIRso(IFILE,recend,tmp, dsize)
C           write(6,*)'rddir1', iout, recend
C
        elseif(recend==recstart) then             
          call rddirso(ifile, recstart, tmp, dsize)
          call myicopy(zlist, tmp(nwrdf+1), nwrdf2-nwrdf)
          CALL WRDIRso(IFILE,recstart,tmp, dsize)
        endif
        deallocate(TMP)
C       ENDIF
C      write(6,*)'========end of putlst======='
      close(415)
300   RETURN
      END
