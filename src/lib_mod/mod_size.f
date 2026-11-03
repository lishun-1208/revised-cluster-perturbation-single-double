! variable must be assignment if use libsrint:
!    iprcln, iprcwd, ichcsz
      module mod_size 
      ! for scalar relativistic CC
      ! ilnbuf only declare in common block
      integer ilnbuf, iprcln, iprcwd ! for length in clmoioso, getlst and putlst
      ! icrsiz only use in readinp and declare in common block
      integer icrsiz,ICHCSZ     !ICRSIZ is the memory size for cpu, ichcsz=4KiB=4096
      ! igpusize only use in readinp, igpuleft not use in libsrint
      integer igpusize, igpuleft !memory size of gpu and memory left in gpu for dynamic calculation

      ! folow variable not use in libsrint
      integer nvvvv, nvvvo, nvvoo, nvooo, noooo, nvvv ! the largest integrals size for different orbitals 
      integer nifvvaa, nifvvbb                                     !fae aa bb
      integer nifvoaa, nifvobb                                     !fme aa bb
      integer nifooaa, nifoobb                                     !fmi aabb
      integer nioi11, nioi12, nioi13                             !iijkl aaaa abab bbbb
      integer nioi21, nioi22, nioi23, nioi24                     !iijka aaaa abab bbbb baba 
      integer nioi31, nioi32, nioi33, nioi34                     !iabij aaaa abab bbbb
      integer nioi41, nioi42, nioi43, nioi44, nioi45, nioi46     !iiabj aaaa abab bbbb baba aabb
      integer nioi51, nioi52, nioi53, nioi54, nioi55             !iabci aaaa abab bbbb baba 
      integer nioi61, nioi62, nioi63                             !iabcd aaaa abab bbbb           
      integer nifvv, nifvo, nifoo, nioi1, nioi2                  !fae, fme and fmi, iijkl, iijka
      integer nioi3, nioi4, nioi5, nioi6                 !iabij, iiabj, iabci, iabcd
      integer icpusize, icpuleft !memory size of gpu and memory left in gpu for dynamic calculation
      integer nioiri_vv, nioiri_vv1, nioiri_vv2
      integer nioiri_vo, nioiri_vo1, nioiri_vo2
      integer nioiri_oo, nioiri_oo1, nioiri_oo2
      end module
        
