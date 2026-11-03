
      INTEGER*8 FUNCTION irpdso(sym,popp,popq,ioppq)
C
C THIS ROUTINE RETURNS THE TOTAL NUMBER OF PQ with symmetry.
C ioppq=0 means p<q and ioppq=1 means p,q
C ioppq=2 means p<=q 
C ioppq = -1, means only one index popp, with popq = 1
C
CEND
      use mod_symm
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION POPP(8), POPQ(8)
      if(sym>8) then 
        write(6,*)'error in irpdso, nirrep>8' 
        stop
      endif
C
      irpdso=0
C      write(6,*)'begin irpdpdso'
      if(ioppq==1) then 
      DO 10 IRREPR=1,NIRREP 
       IRREPL=DIRPRD(IRREPR,sym)
       irpdso=irpdso+popp(irrepl)*popq(irrepr)
10    CONTINUE
      elseif(ioppq==0) then 
      DO 20 IRREPR = 1, NIRREP 
       IRREPL = DIRPRD(IRREPR, SYM) 
       if(irrepl>irrepr) cycle
       if(irrepl<irrepr) then 
         irpdso=irpdso+popp(irrepl)*popq(irrepr)
       elseif(irrepl==irrepr) then 
         irpdso=irpdso+popp(irrepl)*(popq(irrepr)-1)/2
       endif
20    CONTINUE
      elseif(ioppq==-1) then 
      irpdso=popp(sym)
      elseif(ioppq==2) then 
      DO 30 IRREPR = 1, NIRREP
       IRREPL = DIRPRD(IRREPR, SYM)
       if(irrepl>irrepr) cycle
       if(irrepl<irrepr) then
         irpdso=irpdso+popp(irrepl)*popq(irrepr)
       elseif(irrepl==irrepr) then
         irpdso=irpdso+popp(irrepl)*(popq(irrepr)+1)/2
       endif
30    CONTINUE
C
      endif
C
      RETURN
      END
