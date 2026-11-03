
      ! =======================================================================
      ! 预计算单个张量的 DPD 大块基址与块内相对偏移表
      ! =======================================================================
      subroutine build_dpd_tables(pop, symx, inumlft, inumrht, iabsoff, iposlft, iposrht)
      use mod_symm
      implicit none
      integer, intent(in)  :: pop(4, 8)            
      integer, intent(in)  :: symx                 
      integer, intent(out) :: inumlft(33), inumrht(33) 
      integer, intent(out) :: iabsoff(33)   
      integer, intent(out) :: iposlft(33,33), iposrht(33,33) 
      
      integer :: irrep, irr_r, irr_l, irrrht, irrlft, irr_l2
      
      inumlft = 0; inumrht = 0
      iabsoff = 0
      iposlft = 0; iposrht = 0

      do irrep = 1, nirrep
         do irr_r = 1, nirrep
            irr_l = dirprd(irr_r, irrep)
            inumlft(irrep) = inumlft(irrep) + pop(1, irr_l) * pop(2, irr_r)
            inumrht(irrep) = inumrht(irrep) + pop(3, irr_l) * pop(4, irr_r)
         end do
      end do

      iabsoff(1) = 0
      do irrrht = 1, nirrep
         irrlft = dirprd(irrrht, symx)
         iabsoff(irrrht+1) = iabsoff(irrrht) + inumlft(irrlft) * inumrht(irrrht)
         
         iposrht(1, irrrht) = 0
         iposlft(1, irrlft) = 0
         
         do irr_r = 1, nirrep
            irr_l = dirprd(irr_r, irrrht)
            iposrht(irr_r+1, irrrht) = iposrht(irr_r, irrrht) + pop(3, irr_l) * pop(4, irr_r) * inumlft(irrlft)
            
            irr_l2 = dirprd(irr_r, irrlft)
            iposlft(irr_r+1, irrlft) = iposlft(irr_r, irrlft) + pop(1, irr_l2) * pop(2, irr_r)
         end do
      end do
      return
      end subroutine build_dpd_tables
