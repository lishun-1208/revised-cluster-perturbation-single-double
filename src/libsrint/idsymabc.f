
      INTEGER FUNCTION IDSYMABC(SYMABC, POPIN) 
C
C THIS ROUTINE RETURNS THE start point of SIZE OF V(ABC) with special
C IRREP  ! just for a<b<c
C
CEND
      use mod_symm
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION POPP(8), POPQ(8), POPIN(8)
      integer, external :: irpdso
C
!     IF(DIRPRD(IINPA,IINPB)/=DIRPRD(IINPC, SYMABC)) THEN 
!        write(6,*)'error of symmetry in isymoffabc'
!     ENDIF
C
      IDSYMABC = 0
      DO IRREPC = 1, NIRREP
         IRREPAB = DIRPRD(SYMABC, IRREPC)
         DO IRREPB = 1, NIRREP 
            IRREPA = DIRPRD(IRREPB, IRREPAB) 
            NUMA = POPIN(IRREPA)
            NUMB = POPIN(IRREPB)
            NUMC = POPIN(IRREPC)
            if(irrepa>irrepb.or.irrepb>irrepc.or.irrepa>irrepc) cycle 
            if(irrepa==irrepb.and.irrepb==irrepc) then 
               IDSYMABC = IDSYMABC + NUMC*(NUMC-1)*(NUMC-2)/6
!          write(6,'(A,i4,i4,i4,i4)')'idsymabc', idsymabc, irrepa, irrepb, irrepc
            elseif(irrepa<irrepb.and.irrepb==irrepc) then
               IDSYMABC = IDSYMABC + NUMA*NUMB*(NUMB-1)/2
!          write(6,'(A,i4,i4,i4,i4)')'idsymabc', idsymabc, irrepa, irrepb, irrepc
            elseif(irrepa==irrepb.and.irrepb<irrepc) then 
               IDSYMABC = IDSYMABC + NUMA*(NUMA-1)/2*NUMC
!          write(6,'(A,i4,i4,i4,i4)')'idsymabc', idsymabc, irrepa, irrepb, irrepc
            elseif(irrepa<irrepb.and.irrepb<irrepc) then 
               IDSYMABC = IDSYMABC + NUMA*NUMB*NUMC
!          write(6,'(A,i4,i4,i4,i4)')'numa    ', numa, numb, numc
!          write(6,'(A,i4,i4,i4,i4)')'idsymabc', idsymabc, irrepa, irrepb, irrepc
            endif
          ENDDO
        ENDDO
C
      RETURN
      END
