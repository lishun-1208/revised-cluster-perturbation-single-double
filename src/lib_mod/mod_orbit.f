      module mod_orbit
      ! attention! new parameters should be added on the last line
      ! do not add new parameters at the end of existing line

      integer POP(8,2) ! occ orbital for symmetry and spin
      integer VRT(8,2) ! vrt orbital for symmetry and spin
      integer VRTA(8), VRTB(8), POPA(8), POPB(8) ! virtual and occupied orbital for alpha and beta orbital
      integer NV(8), NO(8) ! total virtual and occupied orbitals for different symmetry 
      !integer :: ZORB(8) ! ZORB(1)=1, zorb(2:8)=0 
      integer, parameter, dimension(8) :: zorb = (/1,0,0,0,0,0,0,0/) ! ZORB(1)=1, zorb(2:8)=0 
      
      ! follow variable not use in libsrint
      integer NOCCO(2), NVRTO(2) ! number of occupied and virtual orbital for spin alpha and beta
      integer NVA, NOA, NVB, NOB ! total number of orbital
      integer ZORB2(8) ! zorb2(1:8)=0
      end module
        
