      SUBROUTINE CLMOIOTMPL(IRREPX, ICORE, 
     &           POPP, POPQ, POPR, POPS, filename, iop)
C-this is for large element.
C-index of the file, 412.
C-IRREPX: the symmetry of the tmpfile W(P,Q,R,S) 
C-ICORE: Memory that need to allocate the space 
C-POPP, Q, R, S: the orbital type of the four index. Like POP and VRT in
C--------------- common block
C-filename: the filename to save the W(P,Q,R,S), one should not that the
Clength of 
C--------------- the filename must be (length=8)
C-IOP: to define the type of Wpqrs. 
C-iop==0, then --for P=Q and R=S like Wmnij and Wabef
C----------------------------R-----------I
C IPQRSAAAA       1      1+8
C IPQRSABAB       2      2+8
C IPQRSAAAB       3      3+8
C IPQRSABAA       4      4+8
C IPQRSAABB       5      5+8
C-iop==1, then ---for P=Q and R/=S like Wmnie and Wabei
C IPQRSAABA       6      6+8
C-iop==2, then ---for P/=Q and R/=S like Wmbej
C IPQRSABBA       7      7+8
C IPQRSBAAA       8      8+8
C
C-iop negative for scalcar relativistic calculation
C-iop== -1,    ---for P=Q and R=S
C-iop== -2,    ---for P=Q and R/=S
C-iop== -3,    ---for P/=Q and R/=S
C IPQRSAAAA       1      
C IPQRSABAB       2  
C IPQRSABBA       3   
C
C-iop =-10, just for assym of I<Ab||Cd> and T2(Ab,Ij) 
C   A<=B, and I<=J

      IMPLICIT NONE
C-------------------
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      INTEGER :: SYMRX, SYMRY, SYMRZ
      COMMON/SYMRXYZ/ SYMRX, SYMRY, SYMRZ
C
      INTEGER :: POP, VRT, NT, NF1, NF2
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      INTEGER :: NSTART, NIRREP, IRREPS1, DIRPRD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS1(255,2),DIRPRD(8,8)
C
      INTEGER :: IRPDPD, ISYTYP, NTOT
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      INTEGER :: ISYMOFF
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C----------------------------------

      integer, intent(in) :: popp(8), popq(8)
      integer, intent(in) :: popr(8), pops(8)
      integer, intent(in) :: iop, irrepx
      double precision, intent(in) :: icore(*)
      character*8, intent(in) :: filename
      LOGICAL DOALL,ABIJ,NOABCD
C
C      INTEGER IRREPR(8), IRREPI(8)
      integer ioppq(8), ioprs(8), i, j
      character*2 tmpc
      character*10 filename2
      character*18 pathname
      character*40 opt
      integer ilnbuf, iprcln, iprcwd
      COMMON /FILSPC/ ILNBUF,IPRCLN,IPRCWD
      integer ilength
      integer irreppq, irreprs, irrepp, irrepq, irrepr, irreps
      integer numpq, numrs, numr, nums, nreads
      integer*8 vstart, s
      integer*8 nwords, nfull, npart, totrec, ifile, irecord 
      integer*8 nsize
      integer*8, external :: idsymsoc, irpdso, isymoffso
C
      write(6,*)'begin large clmoiotmp', iprcln
C   
C      call creatfolder(filename)
C      OPEN(UNIT=412, FILE=filename, FORM='UNFORMATTED', 
C     &            access='direct', recl=iprcln)
C      rewind 412
C      write(6,*)'begin izero'
      opt(1:9)='mkdir -p ' 
      opt(10:17)=filename
      write(6,*)'test', trim(opt)
      call system(adjustl(trim(opt(1:17))))  ! to create folderr      
      pathname(1:2)='./'
      pathname(3:10)=filename
      pathname(11:11)='/' 
      CALL IZERO(MOIONUML, 200)
      CALL IZERO(MOIONUMR, 200)
      CALL IZERO(rectmp, 200)
      write(6,*)'clmoioso by wz'
      ifile=1000
C***********************************************************************
      if(iop==1) then
      pathname(12:12)='1'
      do irreprs = 1, nirrep 
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, 0) 
         numrs = irpdso(irreprs, popr, pops, 0) 
         write(pathname(13:13), '(i1)') irreprs
         do irreps = 1, nirrep 
            write(pathname(14:14), '(i1)') irreps
            irrepr = dirprd(irreps, irreprs) 
            numr = popr(irrepr) 
            nums = pops(irreps)
            vstart = isymoffso(irreps, irreprs, popr, pops, 0) 
            if(irrepr>irreps) cycle 
            if(irrepr<irreps) then 
              do s = 1, nums
                 nreads = numr
                 nwords = numr*numpq
                 nfull = nwords/iprcwd
                 npart = mod(nwords, iprcwd)
                 totrec = 1 + nfull
                 write(pathname(15:18), '(i4.4)') s
                 ifile = ifile+1
                 write(6,*)'irrep', ifile, totrec
                 open(unit=ifile, file=pathname, 
     &               form='unformatted',access='direct', recl=iprcln) 
                 do irecord = 1, totrec 
                    call wrdirtmp(ifile, irecord, icore, iprcwd)
                 enddo
                 close(ifile)
              enddo   
            elseif(irrepr==irreps) then 
              do s = 1, nums
                 nreads = s-1
                 nwords = nreads*numpq
                 nfull = nwords/iprcwd
                 npart = mod(nwords, iprcwd)
                 totrec = 1 + nfull
                 write(pathname(15:18), '(i4.4)') s
                 ifile = ifile + 1 
                 write(6,*)'ifile',  ifile, totrec
                 open(unit=ifile, file=trim(pathname), 
     &               form='unformatted',access='direct', recl=iprcln)
                 do irecord = 1, totrec
C                 write(6,*)'ifile', ifile
                    call wrdirtmp(ifile, irecord, icore, iprcwd)
                 enddo
                 close(ifile) 
              enddo
            endif 
         enddo
      enddo
C-----------------------------------------
      pathname(12:12)='2'
      do irreprs = 1, nirrep
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, 0)
         numrs = irpdso(irreprs, popr, pops, 0)
         write(pathname(13:13), '(i1)') irreprs
         do irreps = 1, nirrep
            write(pathname(14:14), '(i1)') irreps
            irrepr = dirprd(irreps, irreprs)
            numr = popr(irrepr)
            nums = pops(irreps)
            vstart = isymoffso(irreps, irreprs, popr, pops, 0) 
              do s = 1, nums
                 nreads = numr
                 nwords = numr*numpq
                 nfull = nwords/iprcwd
                 npart = mod(nwords, iprcwd)
                 totrec = 1 + nfull
                 write(pathname(15:18), '(i4.4)') s
                 ifile = ifile+1
                 write(6,*)'ifile', ifile, totrec 
                 open(unit=ifile, file=pathname, 
     &               form='unformatted',access='direct', recl=iprcln) 
                 do irecord = 1, totrec 
                    call wrdirtmp(ifile, irecord, icore, iprcwd)
                 enddo
                 close(ifile)
              enddo
         enddo
      enddo
C***********************************************************************
      elseif(iop==2) then
      pathname(12:12)='1'
      do irreprs = 1, nirrep
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, 0)
         numrs = irpdso(irreprs, popr, pops, 0)
         write(pathname(13:13), '(i1)') irreprs
         do irreps = 1, nirrep
            write(pathname(14:14), '(i1)') irreps
            irrepr = dirprd(irreps, irreprs)
            numr = popr(irrepr)
            nums = pops(irreps)
            vstart = isymoffso(irreps, irreprs, popr, pops, 0)
            if(irrepr>irreps) cycle
            if(irrepr<irreps) then
              do s = 1, nums
                 nreads = numr
                 nwords = numr*numpq
                 nfull = nwords/iprcwd
                 npart = mod(nwords, iprcwd)
                 totrec = 1 + nfull
                 write(pathname(15:18), '(i4.4)') s
                 ifile = ifile+1
                 write(6,*)'irrep', ifile, totrec
                 open(unit=ifile, file=pathname,
     &               form='unformatted',access='direct', recl=iprcln)
                 do irecord = 1, totrec
                    call wrdirtmp(ifile, irecord, icore, iprcwd)
                 enddo
                 close(ifile)
              enddo
            elseif(irrepr==irreps) then
              do s = 1, nums
                 nreads = s
                 nwords = nreads*numpq
                 nfull = nwords/iprcwd
                 npart = mod(nwords, iprcwd)
                 totrec = 1 + nfull
                 write(pathname(15:18), '(i4.4)') s
                 ifile = ifile + 1
                 write(6,*)'ifile',  ifile, totrec
                 open(unit=ifile, file=trim(pathname),
     &               form='unformatted',access='direct', recl=iprcln)
                 do irecord = 1, totrec
C                 write(6,*)'ifile', ifile
                    call wrdirtmp(ifile, irecord, icore, iprcwd)
                 enddo
                 close(ifile)
              enddo
            endif
         enddo
      enddo
      endif
C-----------------------------------------
C      OPEN(UNIT=413, FILE=filename2,
C     &                 FORM='UNFORMATTED',status='unknown')
C      rewind 413 
C      write(413) IRREPX
C      write(413) (popp(i), i=1,8)
C      write(413) (popq(i), i=1,8)
C      write(413) (popr(i), i=1,8)
C     write(413) (pops(i), i=1,8)
C      write(413) (ioppq(i), i = 1,8)
C      write(413) (ioprs(i), i = 1,8)
C      write(413) (irrepr(i), i = 1,8)
C      write(413) (irrepi(i), i = 1,8)
C      write(413) ((moionuml(i,j),i=1,10),j=1,20)
C      write(413) ((moionumr(i,j),i=1,10),j=1,20)
C      write(413) ((rectmp(i,j),i=1,10),j=1,20)
C      write(413) iop
C      close(413)
C      write(6,*)'end of file',' ',filename2
C      stop
C***********************************************************************
      RETURN
      END
