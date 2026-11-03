      SUBROUTINE PUTLSTTMPL(ZLIST,idx, irreps, irreprs, ilist,
     &                           LENGTH,filename)
C
CEND
      IMPLICIT NONE
      integer, intent(in)::ilist, irreprs, irreps, idx, length 
      character(len=8), intent(in) :: filename
C
      double precision, intent(inout) :: ZLIST(*) 
C     
C      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
C     &               MOIODS(10,500),MOIOFL(10,500)
      INTEGER ILNBUF, IPRCLN, IPRCWD
      COMMON /FILSPC/ ILNBUF,IPRCLN,IPRCWD
C
      INTEGER*8 ICRSIZ, ICHCSZ, IOFFX, LENREC
      COMMON /IOPOS/ ICRSIZ,ICHCSZ,IOFFX(2),LENREC
C
      INTEGER IFLAGS
      COMMON /FLAGS/ IFLAGS(100)
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      double precision, allocatable :: TMP(:)
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
      character(len=18) pathname
C
      integer*8 nwords, nstart, nend, nwrdf, nwrdf2
      integer nfull, npart
      integer irrepx, i, j, ifile, listrec
      integer recstart, recend, istat
      integer iout, irecord, istart, totrec
      integer*8 tmpi1, tmpi2
C
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
C      goto 300
      nwords = length 
      nfull = nwords/iprcwd
      npart = mod(nwords, iprcwd) 
      totrec = 1 + nfull 
C
      pathname(1:2)='./'
      pathname(3:10)=filename
      pathname(11:11)='/'
      write(pathname(12:12), '(i1)') ilist
      write(pathname(13:13), '(i1)') irreprs 
      write(pathname(14:14), '(i1)') irreps
      write(pathname(15:18), '(i4.4)') idx 
     
      open(unit=412, file=trim(pathname),
     &               form='unformatted',access='direct', recl=iprcln)
      do irecord = 1, totrec-1
         istart = (irecord-1)*iprcwd+1 
         call wrdirso(412, irecord, zlist(istart), iprcwd)
      enddo
         istart=(totrec-1)*iprcwd+1
         call wrdirso(412, totrec, zlist(istart), iprcwd)
      close(412)
C
C DEAL WITH RIDICULOUS LIST NUMBERS
C
C      write(6,*)'========end of putlst======='
300   RETURN
      END
