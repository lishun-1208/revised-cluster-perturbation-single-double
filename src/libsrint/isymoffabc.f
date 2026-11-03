
      INTEGER FUNCTION ISYMOFFABC(IINPA, IINPB, IINPC, IRREPABC, POPIN) 
C
C THIS ROUTINE RETURNS THE start point of SIZE OF V(ABC) with special
C IRREP  ! just fo a<b<c
C
CEND
      use mod_symm
      IMPLICIT NONE 
C
      integer iinpa, iinpb, iinpc, irrepabc, IENDB
      INTEGER IRREPC, IRREPB, IRREPA, irrepab
      INTEGER POPP(8), POPQ(8), POPIN(8)
      integer numa, numb, numc
      integer, external :: irpdso
C
      IF(DIRPRD(IINPA,IINPB)/=DIRPRD(IINPC, irrepabc)) THEN 
         write(6,*)'error of symmetry in isymoffabc'
         write(6,*) iinpa, iinpb, iinpc, irrepabc 
      ENDIF
C
C     IF(iinpa>iinpb.or.iinpb>iinpc.or.iinpa>iinpc) THEN 
C        write(6,*)'error of iinpabc in isymoffabc'
C     ENDIF
      ISYMOFFABC = 1
      DO IRREPC = 1, IINPC
         IRREPAB = DIRPRD(IRREPABC, IRREPC)
         IENDB = NIRREP
         IF(IRREPC==IINPC) IENDB = IINPB-1
         DO IRREPB = 1, IENDB
            IRREPA = DIRPRD(IRREPB, IRREPAB) 
            NUMA = POPIN(IRREPA)
            NUMB = POPIN(IRREPB)
            NUMC = POPIN(IRREPC)
            if(irrepa>irrepb.or.irrepb>irrepc.or.irrepa>irrepc) cycle 
            if(irrepa==irrepb.and.irrepb==irrepc) then 
               ISYMOFFABC = ISYMOFFABC + NUMC*(NUMC-1)*(NUMC-2)/6
            elseif(irrepa<irrepb.and.irrepb==irrepc) then
               ISYMOFFABC = ISYMOFFABC + NUMA*NUMB*(NUMB-1)/2
            elseif(irrepa==irrepb.and.irrepb<irrepc) then 
               ISYMOFFABC = ISYMOFFABC + NUMA*(NUMA-1)/2*NUMC
            elseif(irrepa<irrepb.and.irrepb<irrepc) then 
               ISYMOFFABC = ISYMOFFABC + NUMA*NUMB*NUMC
            endif
          ENDDO
        ENDDO
C
      RETURN
      END
