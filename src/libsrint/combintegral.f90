 subroutine combintegral(iout, iin, poppa, poppb, popqa, popqb, ioppq, & 
                          popra, poprb, popsa, popsb, ioprs, symout, & 
                          poppn, popqn, ioppqn, poprn, popsn, ioprsn, symin, intype, itest)
! combine spin integral to nspin integral.
! iin for spin integral, iout for nspin integral
! intype for close shell
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA
!
! spin case need to be put in nspin integral time-reversal symmetry must be employed
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB! ABBB, BBAA, BBAB, BBBB 
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, ABBB, BBAA, BBAB, BBBB, AABA, ABBA, BBBA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA, BBAA, BBAB, BBBB, BAAA, BAAB, BABB
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA, BBBB, 
!                    AABA, ABBA, BBBA, BAAA, BAAB, BABB, BABA
! example::
! example:: T2(abij) ABIJ, ABIj, AbIJ, ABij, AbIj
! example:: R2(abci) ABCI, ABCi, AbCI, ABci, AbCi, ABcI
! example:: W2(mbej) MBEJ, MBEj, MbEJ, MBej, MbEj, MBeJ, mBEJ, MbeJ
!     ioppq_spin = 0; ioprs_spin = 1; sym_spin = symry
!     ioppq_ns   = 0; ioprs_ns   = 0; sym_nspin = 1 
!     call zeroso(t2r,idsymsoc(sym_spin,nv,nv,ioppq_ns,no,no,ioprs_ns)
!     call combintegral(t2r, t2aaabr, vrta, vrtb, vrta, vrtb, ioppq_spin, &
!                        popa, popb, popa, popb, ioprs_spin, sym_spin, &
!                        nv, nv, 0, no, no, 0, sym_nspin, 'AAAB')
! type for open shell
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, ABBB, BBAA, BBAB, BBBB
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, ABBB, BBAA, BBAB, BBBB, AABA, ABBA, BBBA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA, BBAA, BBAB, BBBB, BAAA, BAAB, BABB
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA, BBBB, 
!                    AABA, ABBA, BBBA, BAAA, BAAB, BABB, BABA
 use mod_iop
 use mod_size
 use mod_orbit
 use mod_ioff
 use mod_symm

 implicit none
!
 integer, intent(in) :: symin, symout, ioppqn, ioprsn, ioppq, ioprs
 integer, intent(in) :: poppa(8), poppb(8), popqa(8), popqb(8) 
 integer, intent(in) :: popra(8), poprb(8), popsa(8), popsb(8)
 integer, intent(in) :: poppn(8), popqn(8), poprn(8), popsn(8)
!
 double precision, intent(inout) :: iin(*)
 double precision, intent(inout) :: iout(*), itest(*)
!
 character*4, intent(inout):: intype
 integer :: popp_s(8), popq_s(8), popr_s(8), pops_s(8)
 integer :: irreppq, irreprs, irrepp, irrepq, irrepr, irreps
 integer :: irreppqn, irreprsn, irreppn, irrepqn, irreprn, irrepsn
 integer :: numpn, numqn, numrn, numsn, numpqn, numrsn
 integer :: p, q, r, s, numpq, numrs, nump, numq, numr, nums
 integer :: ioffns(8), ioffnspin, ioffspin
 integer :: irrepqp, irrepsr, numqp, numsr 
 integer :: ioffin(9), ioffout(9), ioffin2(9)
 integer :: ioffp, ioffq, ioffr, ioffs, idx, idxn
 integer :: nsize, sym2, idxn2
 integer :: numpqn2, numrsn2, ioffns2(8), ioffnspin2
 character*4 str
!
 integer, external :: irpdso, idsymsoc, isymoffso
 real*8, external :: nrm2so 
!
 if(prec==8) then 
   call combintegral_d(iout, iin, poppa, poppb, popqa, popqb, ioppq, & 
                          popra, poprb, popsa, popsb, ioprs, symout, & 
                          poppn, popqn, ioppqn, poprn, popsn, ioprsn,&
                          symin, intype, itest)
 elseif(prec==4) then 
   call combintegral_s(iout, iin, poppa, poppb, popqa, popqb, ioppq, & 
                          popra, poprb, popsa, popsb, ioprs, symout, & 
                          poppn, popqn, ioppqn, poprn, popsn, ioprsn,&
                          symin, intype, itest)
 endif
 return
 ioffout(1) = 1; ioffin(1) = 1
 ioffns(1) = 1; ioffns2(1) = 1
 ioffin2(1) = 1
 do irreprs = 1, nirrep 
    irreppq = dirprd(irreprs, symout)
     if(intype=='AAAA') then 
!--aaaa
       numpq = irpdso(irreppq, poppa, popqa, ioppq) 
       numrs = irpdso(irreprs, popra, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqa; popr_s = popra; pops_s = popsa
 !      write(6,*)'pop', popp_s, poppa
    elseif(intype=='ABAB') then  
!--abab
       numpq = irpdso(irreppq, poppa, popqb, ioppq)
       numrs = irpdso(irreprs, popra, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqb; popr_s = popra; pops_s = popsb
    elseif(intype=='AABB') then  
!--aabb
       numpq = irpdso(irreppq, poppa, popqa, ioppq) 
       numrs = irpdso(irreprs, poprb, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqa; popr_s = poprb; pops_s = popsb
    elseif(intype=='AAAB') then  
!--aaab
       numpq = irpdso(irreppq, poppa, popqa, ioppq) 
       numrs = irpdso(irreprs, popra, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqa; popr_s = popra; pops_s = popsb
    elseif(intype=='ABAA') then  
!--abaa 
       numpq = irpdso(irreppq, poppa, popqb, ioppq) 
       numrs = irpdso(irreprs, popra, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqb; popr_s = popra; pops_s = popsa
!==========end close shell for P=Q and R=S
    elseif(intype=='ABBB') then  
!--abbb 
       numpq = irpdso(irreppq, poppa, popqb, ioppq) 
       numrs = irpdso(irreprs, poprb, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqb; popr_s = poprb; pops_s = popsb
    elseif(intype=='BBAA') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqb, ioppq) 
       numrs = irpdso(irreprs, popra, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqb; popr_s = popra; pops_s = popsa
    elseif(intype=='BBAB') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqb, ioppq) 
       numrs = irpdso(irreprs, popra, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqb; popr_s = popra; pops_s = popsb
    elseif(intype=='BBBB') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqb, ioppq) 
       numrs = irpdso(irreprs, poprb, popsb, ioprs)
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqb; popr_s = poprb; pops_s = popsb
!==========end open shell for P=Q and R=S
    elseif(intype=='AABA') then  
!--abbb 
       numpq = irpdso(irreppq, poppa, popqa, ioppq) 
       numrs = irpdso(irreprs, poprb, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqa; popr_s = poprb; pops_s = popsa
    elseif(intype=='ABBA') then  
!--abbb 
       numpq = irpdso(irreppq, poppa, popqb, ioppq) 
       numrs = irpdso(irreprs, poprb, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppa; popq_s = popqb; popr_s = poprb; pops_s = popsa
    elseif(intype=='BBBA') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqb, ioppq) 
       numrs = irpdso(irreprs, poprb, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqb; popr_s = poprb; pops_s = popsa
!========== for P=Q
    elseif(intype=='BAAA') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqa, ioppq) 
       numrs = irpdso(irreprs, popra, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqa; popr_s = popra; pops_s = popsa
    elseif(intype=='BAAB') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqa, ioppq) 
       numrs = irpdso(irreprs, popra, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqa; popr_s = popra; pops_s = popsb
    elseif(intype=='BABB') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqa, ioppq) 
       numrs = irpdso(irreprs, poprb, popsb, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqa; popr_s = poprb; pops_s = popsb
!========== for R=S 
    elseif(intype=='BABA') then  
!--abbb 
       numpq = irpdso(irreppq, poppb, popqa, ioppq) 
       numrs = irpdso(irreprs, poprb, popsa, ioprs) 
       ioffout(irreprs+1) = ioffout(irreprs) + numpq*numrs
       popp_s = poppb; popq_s = popqa; popr_s = poprb; pops_s = popsa
!========== for independent P Q R S
    else
       write(6,*)'error in split3 ioff'
       stop
    endif
!--ns
    irreppq = dirprd(irreprs, symin)
    numpq = irpdso(irreppq, poppn, popqn, ioppqn) 
    numrs = irpdso(irreprs, poprn, popsn, ioprsn)
    ioffin(irreprs+1) = ioffin(irreprs) + numpq*numrs
    numpq = irpdso(irreppq, poppn, popqn, 1) 
    numrs = irpdso(irreprs, poprn, popsn, 1)
    ioffin2(irreprs+1) = ioffin2(irreprs) + numpq*numrs
 enddo            
   nsize = idsymsoc(1, poppn, popqn, 1, poprn, popsn, 1)
!-aaaa
   do irreprs = 1, nirrep 
      irreppq = dirprd(irreprs, symout) 
      numrs = irpdso(irreprs, popr_s, pops_s, ioprs) 
      numpq = irpdso(irreppq, popp_s, popq_s, ioppq)
!      write(6,*)'num', numrs, numpq 
      do irreps = 1, nirrep 
         irrepr = dirprd(irreps, irreprs) 
         nums = pops_s(irreps) 
         numr = popr_s(irrepr) 
         do irrepq = 1, nirrep 
            irrepp = dirprd(irrepq, irreppq) 
            numq = popq_s(irrepq) 
            nump = popp_s(irrepp)
            ioffspin = ioffout(irreprs) -1 & 
                + (isymoffso(irreps, irreprs, popr_s, pops_s, ioprs)-1)*numpq & 
                +  isymoffso(irrepq, irreppq, popp_s, popq_s, ioppq) -1
                if(intype=='AAAA') then 
                  irreppn = irrepp
                  irrepqn = irrepq
                  irreprn = irrepr
                  irrepsn = irreps
                  ioffp = 0
                  ioffq = 0
                  ioffr = 0
                  ioffs = 0               
                elseif(intype=='ABAB') then 
                  irreppn = irrepp 
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = irrepr
                  irrepsn = dirprd(irreps, symry)
                  ioffp = 0
                  ioffq = popqa(irrepqn) 
                  ioffr = 0
                  ioffs = popsa(irrepsn)
                elseif(intype=='AABB') then 
                  irreppn = irrepp 
                  irrepqn = irrepq
                  irreprn = dirprd(irrepr, symry) 
                  irrepsn = dirprd(irreps, symry)
                  ioffp = 0
                  ioffq = 0
                  ioffr = popra(irreprn) 
                  ioffs = popsa(irrepsn)
                elseif(intype=='AAAB') then 
                  irreppn = irrepp 
                  irrepqn = irrepq 
                  irreprn = irrepr
                  irrepsn = dirprd(irreps, symry)
                  ioffp = 0
                  ioffq = 0 
                  ioffr = 0
                  ioffs = popsa(irrepsn)
                elseif(intype=='ABAA') then 
                  irreppn = irrepp 
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = irrepr
                  irrepsn = irreps                 
                  ioffp = 0
                  ioffq = popqa(irrepqn) 
                  ioffr = 0
                  ioffs = 0  
                elseif(intype=='ABBB') then 
                  irreppn = irrepp 
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = dirprd(irrepr, symry)
                  irrepsn = dirprd(irreps, symry)               
                  ioffp = 0
                  ioffq = popqa(irrepqn) 
                  ioffr = popra(irreprn)
                  ioffs = popsa(irrepsn)
                elseif(intype=='ABBA') then 
                  irreppn = irrepp 
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = dirprd(irrepr, symry)
                  irrepsn = irreps               
                  ioffp = 0
                  ioffq = popqa(irrepqn) 
                  ioffr = popra(irreprn)
                  ioffs = 0 
                elseif(intype=='BBAA') then 
                  irreppn = dirprd(irrepp, symry)
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = irrepr
                  irrepsn = irreps            
                  ioffp = poppa(irreppn)
                  ioffq = popqa(irrepqn) 
                  ioffr = 0 
                  ioffs = 0
                elseif(intype=='BBAB') then 
                  irreppn = dirprd(irrepp, symry)
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = irrepr
                  irrepsn = dirprd(irreps, symry)
                  ioffp = poppa(irreppn)
                  ioffq = popqa(irrepqn) 
                  ioffr = 0 
                  ioffs = popsa(irrepsn)
                elseif(intype=='BBBB') then 
                  irreppn = dirprd(irrepp, symry)
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = dirprd(irrepr, symry)
                  irrepsn = dirprd(irreps, symry)
                  ioffp = poppa(irreppn)
                  ioffq = popqa(irrepqn) 
                  ioffr = popra(irreprn)
                  ioffs = popsa(irrepsn)
                elseif(intype=='BBBA') then 
                  irreppn = dirprd(irrepp, symry)
                  irrepqn = dirprd(irrepq, symry) 
                  irreprn = dirprd(irrepr, symry)
                  irrepsn = irreps            
                  ioffp = poppa(irreppn)
                  ioffq = popqa(irrepqn) 
                  ioffr = popra(irreprn)
                  ioffs = 0
                elseif(intype=='BAAA') then 
                  irreppn = dirprd(irrepp,symry)
                  irrepqn = irrepq
                  irreprn = irrepr
                  irrepsn = irreps
                  ioffp = poppa(irreppn)
                  ioffq = 0
                  ioffr = 0
                  ioffs = 0               
                elseif(intype=='BAAB') then 
                  irreppn = dirprd(irrepp,symry)
                  irrepqn = irrepq
                  irreprn = irrepr
                  irrepsn = dirprd(irreps,symry)
                  ioffp = poppa(irreppn)
                  ioffq = 0
                  ioffr = 0
                  ioffs = popsa(irrepsn)  
                elseif(intype=='BABB') then 
                  irreppn = dirprd(irrepp,symry)
                  irrepqn = irrepq
                  irreprn = dirprd(irrepr,symry)
                  irrepsn = dirprd(irreps,symry)
                  ioffp = poppa(irreppn)
                  ioffq = 0
                  ioffr = popra(irreprn)
                  ioffs = popsa(irrepsn)  
                elseif(intype=='BABA') then 
                  irreppn = dirprd(irrepp,symry)
                  irrepqn = irrepq
                  irreprn = dirprd(irrepr,symry)
                  irrepsn = irreps              
                  ioffp = poppa(irreppn)
                  ioffq = 0
                  ioffr = popra(irreprn)
                  ioffs = 0               
                elseif(intype=='AABA') then 
                  irreppn = irrepp 
                  irrepqn = irrepq 
                  irreprn = dirprd(irrepr, symry)
                  irrepsn = irreps                 
                  ioffp = 0
                  ioffq = 0 
                  ioffr = popra(irreprn)
                  ioffs = 0  
                else 
                  write(6,*)'error in split3', intype
                  stop           
                endif
            irreprsn = dirprd(irreprn, irrepsn) 
            irreppqn = dirprd(irreppn, irrepqn) 
            numpqn = irpdso(irreppqn, poppn, popqn, ioppqn) 
            numrsn = irpdso(irreprsn, poprn, popsn, ioprsn)
            numpn = poppn(irreppn)
            numqn = popqn(irrepqn) 
            numrn = poprn(irreprn) 
            numsn = popsn(irrepsn)
!          if(intype=='AAAA') then
!             !write(6,*)'iop', ioppq, ioprs, ioppqn, ioprsn
!             call splitintegral_sub3(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, & 
!                   irrepp, irrepq, irrepr, irreps, ioffspin, iout, &
!                   poppn, popqn, ioppqn, poprn, popsn, ioprsn, &
!                   ioffp, ioffq, ioffr, ioffs, &
!                   irreppn, irrepqn, irreprn, irrepsn, ioffin, iin)
!          elseif(intype=='AAAB') then
!             call splitintegral_sub3(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, & 
!                   irrepp, irrepq, irrepr, irreps, ioffspin, iout, &
!                   poppn, popqn, ioppqn, poprn, popsn, ioprsn, &
!                   ioffp, ioffq, ioffr, ioffs, &
!                   irreppn, irrepqn, irreprn, irrepsn, ioffin, iin)
!          elseif(intype=='ABAA') then
!             call splitintegral_sub3(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, & 
!                   irrepp, irrepq, irrepr, irreps, ioffspin, iout, &
!                   poppn, popqn, ioppqn, poprn, popsn, ioprsn, &
!                   ioffp, ioffq, ioffr, ioffs, &
!                   irreppn, irrepqn, irreprn, irrepsn, ioffin, iin)
!          elseif(intype=='AABB') then
!             call splitintegral_sub3(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, & 
!                   irrepp, irrepq, irrepr, irreps, ioffspin, iout, &
!                   poppn, popqn, ioppqn, poprn, popsn, ioprsn, &
!                   ioffp, ioffq, ioffr, ioffs, &
!                   irreppn, irrepqn, irreprn, irrepsn, ioffin, iin)
!          elseif(intype=='ABAB') then
!             if(intype=='BBBB') then 
!                write(6,*)'ioffspin', irreprs, ioffout(irreprs)
!             endif
              call combintegral_sub(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, & 
                    irrepp, irrepq, irrepr, irreps, ioffspin, iin, &
                    poppn, popqn, ioppqn, poprn, popsn, ioprsn, &
                    ioffp, ioffq, ioffr, ioffs, &
                    irreppn, irrepqn, irreprn, irrepsn, ioffin, iout)
!          endif
         enddo
      enddo
    enddo
!   if(intype=='AAAB') then
!   write(6,*)'norm of t2aaaar ', nrm2so(idsymsoc(1,poppn,popqn,0,poprn,popsn,0),iin,1)*2.d0
!   write(6,*)'norm of t2aaaar ', nrm2so(idsymsoc(1,poppn,popqn,1,poprn,popsn,1),itest,1)
!   do irreprsn = 1, nirrep 
!      irreppqn = dirprd(irreprsn, 1) 
!      numrsn  = irpdso(irreprsn, poprn, popsn, 0)
!      numpqn  = irpdso(irreppqn, poppn, popqn, 0)
!      numrsn2 = irpdso(irreprsn, poprn, popsn, 1)
!      numpqn2 = irpdso(irreppqn, poppn, popqn, 1)
!      do irrepsn = 1, nirrep
!         irreprn = dirprd(irrepsn, irreprsn) 
!         do irrepqn = 1, nirrep 
!            irreppn = dirprd(irrepqn, irreppqn) 
!            numpn = poppn(irreppn) 
!            numqn = popqn(irrepqn) 
!            numrn = poprn(irreprn) 
!            numsn = popsn(irrepsn) 
!            ioffnspin = ioffin(irreprsn) - 1  &
!              + (isymoffso(irrepsn, irreprsn, poprn, popsn,0)-1)*numpqn &
!              +  isymoffso(irrepqn, irreppqn, poppn, popqn,0) - 1
!            ioffnspin2 = ioffin2(irreprsn) - 1  &
!              + (isymoffso(irrepsn, irreprsn, poprn, popsn,1)-1)*numpqn2 &
!              +  isymoffso(irrepqn, irreppqn, poppn, popqn,1) - 1
!            if(irreppn<irrepqn.and.irreprn<irrepsn) then 
!               do s = 1, numsn
!               do r = 1, numrn
!               do q = 1, numqn
!               do p = 1, numpn
!                  idxn = ioffnspin + ((s-1)*numrn+r-1)*numpqn & 
!                            + (q-1)*numpn + p
!                  idxn2 = ioffnspin2 + ((s-1)*numrn+r-1)*numpqn2 & 
!                            + (q-1)*numpn + p
!       !          write(6,*)'idx', idxn, idxn2
!                  write(6,'(A,2I10,2F18.12)')'iintest2', idxn, idxn2, iin(idxn), itest(idxn2)
!               enddo
!               enddo
!               enddo
!               enddo
!            endif
!         enddo
!      enddo
!   enddo            
!   endif
    return
    end
