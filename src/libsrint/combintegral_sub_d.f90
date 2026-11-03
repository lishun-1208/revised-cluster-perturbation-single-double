 subroutine combintegral_sub_d(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, &
                 irrepp, irrepq, irrepr, irreps, ioffspin, iout,         &        
                 poppn, popqn, ioppqn, poprn, popsn, ioprsn,   &
                 ioffp, ioffq, ioffr, ioffs,                    &
                 irreppn, irrepqn, irreprn, irrepsn, ioffin, iin)      
 use mod_symm
 implicit none
 integer, intent(in):: nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, ioffspin
 integer, intent(in):: poppn(8), popqn(8), ioppqn, poprn(8), popsn(8), ioprsn 
 integer, intent(in):: ioffp, ioffq, ioffr, ioffs,  ioffin(8) 
 integer, intent(in):: irrepp, irrepq, irrepr, irreps
 integer, intent(in) :: irreppn, irrepqn, irreprn, irrepsn
 real*8, intent(in) :: iin(*)
 real*8, intent(out) :: iout(*)
 integer ioppqn_calc, ioprsn_calc, ioppq_calc, ioprs_calc, icpqn, icrsn
 integer ioffnspin, irreppqn, irreprsn
 integer numpn, numqn, numrn, numsn, numpqn, numrsn
  
 integer, external :: irpdso, isymoffso

 numpn = poppn(irreppn)
 numqn = popqn(irrepqn)
 numrn = poprn(irreprn)
 numsn = popsn(irrepsn)
 irreppqn = dirprd(irreppn, irrepqn) 
 irreprsn = dirprd(irreprn, irrepsn) 
 numpqn = irpdso(irreppqn, poppn, popqn, ioppqn) 
 numrsn = irpdso(irreprsn, poprn, popsn, ioprsn)
 if(ioppqn==0.and.ioprsn==0) then 
    if(irreppn<irrepqn.and.irreprn<irrepsn) then
      ioppqn_calc=1; ioprsn_calc=1
      icpqn = 0    ; icrsn      =0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1

    elseif(irreppn<irrepqn.and.irreprn==irrepsn) then
      ioppqn_calc=1; ioprsn_calc=0
      icpqn = 0    ; icrsn      =0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1

    elseif(irreppn<irrepqn.and.irreprn>irrepsn) then
      ioppqn_calc=1; ioprsn_calc=1
      icpqn = 0    ; icrsn      =1
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irreprn, irreprsn, popsn, poprn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1

    elseif(irreppn==irrepqn.and.irreprn<irrepsn) then
      ioppqn_calc=0; ioprsn_calc=1
      icpqn = 0    ; icrsn      =0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1

    elseif(irreppn==irrepqn.and.irreprn==irrepsn) then
      ioppqn_calc=0; ioprsn_calc=0
      icpqn = 0    ; icrsn      =0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1

    elseif(irreppn==irrepqn.and.irreprn>irrepsn) then
      ioppqn_calc=0; ioprsn_calc=1
      icpqn = 0    ; icrsn      =1
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irreprn, irreprsn, popsn, poprn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1

    elseif(irreppn>irrepqn.and.irreprn<irrepsn) then
      ioppqn_calc=1; ioprsn_calc=1
      icpqn = 1    ; icrsn      =0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irreppn, irreppqn, popqn, poppn,ioppqn) - 1

    elseif(irreppn>irrepqn.and.irreprn==irrepsn) then
      ioppqn_calc=1; ioprsn_calc=0
      icpqn = 1    ; icrsn      =0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irreppn, irreppqn, popqn, poppn,ioppqn) - 1

    elseif(irreppn>irrepqn.and.irreprn>irrepsn) then
      ioppqn_calc=1; ioprsn_calc=1
      icpqn = 1    ; icrsn      =1
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irreprn, irreprsn, popsn, poprn,ioprsn)-1)*numpqn&
       +  isymoffso(irreppn, irreppqn, popqn, poppn,ioppqn) - 1
    else
      write(6,*)'error in splitintegral_sub3'
      stop
    endif
 elseif(ioppqn==0.and.ioprsn==1) then
    icrsn = 0 
    if(irreppn<irrepqn) then
      ioppqn_calc=1; ioprsn_calc=1; icpqn=0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1
    elseif(irreppn==irrepqn) then
      ioppqn_calc=0; ioprsn_calc=1; icpqn=0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1
    elseif(irreppn>irrepqn) then
      ioppqn_calc=1; ioprsn_calc=1; icpqn=1
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irreppn, irreppqn, popqn, poppn,ioppqn) - 1
    else 
      return
    endif
 elseif(ioppqn==1.and.ioprsn==0) then
    icpqn = 0 
    if(irreprn<irrepsn) then
      ioppqn_calc=1; ioprsn_calc=1; icrsn=0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1
    elseif(irreprn==irrepsn) then
      ioppqn_calc=1; ioprsn_calc=0; icrsn=0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1
    elseif(irreprn>irrepsn) then
      ioppqn_calc=1; ioprsn_calc=1; icrsn=1
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irreprn, irreprsn, popsn, poprn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1
    else
      return
    endif
 elseif(ioppqn==1.and.ioprsn==1) then 
      ioppqn_calc=1; ioprsn_calc=1
      icpqn = 0; icrsn = 0
      ioffnspin = ioffin(irreprsn) - 1 &
       + (isymoffso(irrepsn, irreprsn, poprn, popsn,ioprsn)-1)*numpqn&
       +  isymoffso(irrepqn, irreppqn, poppn, popqn,ioppqn) - 1
 endif

 if(ioppq==0.and.ioprs==0) then 
    if(irrepp<irrepq.and.irrepr<irreps) then
      ioppq_calc=1; ioprs_calc=1
    elseif(irrepp==irrepq.and.irrepr<irreps) then
      ioppq_calc=0; ioprs_calc=1
    elseif(irrepp<irrepq.and.irrepr==irreps) then
      ioppq_calc=1; ioprs_calc=0
    elseif(irrepp==irrepq.and.irrepr==irreps) then
      ioppq_calc=0; ioprs_calc=0
    else
      return
    endif
 elseif(ioppq==0.and.ioprs==1) then 
    if(irrepp<irrepq) then
      ioppq_calc=1; ioprs_calc=1
    elseif(irrepp==irrepq) then
      ioppq_calc=0; ioprs_calc=1
    else
      return
    endif
 elseif(ioppq==1.and.ioprs==0) then 
    if(irrepr<irreps) then
      ioppq_calc=1; ioprs_calc=1
    elseif(irrepr==irreps) then
      ioppq_calc=1; ioprs_calc=0
    else
      return
    endif
 elseif(ioppq==1.and.ioprs==1) then 
      ioppq_calc=1; ioprs_calc=1
 endif
!write(6,'(A, 5I15)')'in2out3 begin', ioffp, ioffq, ioffr, ioffs, ioffnspin
!write(6,'(A, I15)')'ioff', ioffin
!write(6,'(A, 6I15)')'num', numpn, numqn, numrn, numsn, numpqn, numrsn
! write(6,'(A, 6I5)')'iop', ioppq_calc, ioprs_calc, ioppqn_calc, ioprsn_calc, icpqn, icrsn
 call combintegral_spin2nspin_d(nump, numq, ioppq_calc, numr, nums, ioprs_calc, &
        numpq, numrs, ioffspin, iout, &
        numpn, numqn, ioppqn_calc, icpqn,numrn, numsn, ioprsn_calc,icrsn, numpqn, numrsn, &
        ioffp, ioffq, ioffr, ioffs, ioffnspin, iin)
 return
 end

  subroutine spin2nspin_d(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, ioffspin, iout, &
                   numpn, numqn, ioppqn, icpqn, numrn, numsn, ioprsn, icrsn, numpqn, numrsn, & 
                   ioffp, ioffq, ioffr, ioffs, ioffnspin, iin)
  implicit none 
  integer,intent(in):: nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, ioffspin
  integer,intent(in):: numpn, numqn, ioppqn, numrn, numsn, ioprsn, numpqn, numrsn 
  integer, intent(in):: ioffp, ioffq, ioffr, ioffs, ioffnspin, icpqn, icrsn
  real*8, intent(inout) :: iin(*)
  real*8, intent(inout) :: iout(*)
  integer p, q, r, s, idx, idxn
 ! write(6,*)'iop', ioppqn, ioprsn, ioppq, ioprs
  !1111111111-----------------------------------------------------------
  if(ioppqn==1.and.ioprsn==1) then 
  !1222222222-----------------------------------------------------------
    if(icpqn==0.and.icrsn==0) then 
      if(ioppq==1.and.ioprs==1) then
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx) 
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==1) then 
         if(numq>1) then              
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx) 
         enddo
         enddo
         enddo
         enddo
         endif
      elseif(ioppq==1.and.ioprs==0) then             
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx) 
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==0) then              
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx) 
         enddo
         enddo
         enddo
         enddo
      endif
  !1122222222-----------------------------------------------------------
    elseif(icpqn==0.and.icrsn==1) then 
      if(ioppq==1.and.ioprs==1) then
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==1) then
         if(numq>1) then
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
         endif
      elseif(ioppq==1.and.ioprs==0) then
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==0) then
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
      endif
  !1133333333-----------------------------------------------------------
    elseif(icpqn==1.and.icrsn==0) then 
      if(ioppq==1.and.ioprs==1) then
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==1) then
         if(numq>1) then
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
         endif
      elseif(ioppq==1.and.ioprs==0) then
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==0) then
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
         enddo
         enddo
         enddo
         enddo
      endif
  !1144444444-----------------------------------------------------------
    elseif(icpqn==1.and.icrsn==1) then 
      if(ioppq==1.and.ioprs==1) then
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = iout(idx)
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==1) then
         if(numq>1) then
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = iout(idx)
         enddo
         enddo
         enddo
         enddo
         endif
      elseif(ioppq==1.and.ioprs==0) then
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = iout(idx)
         enddo
         enddo
         enddo
         enddo
      elseif(ioppq==0.and.ioprs==0) then
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = iout(idx)
         enddo
         enddo
         enddo
         enddo
      endif
     endif

  !2222222222-----------------------------------------------------------
  elseif(ioppqn==0.and.ioprsn==1) then
    if(icpqn==0.and.icrsn==0) then  
     if(ioppq==1.and.ioprs==1) then
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==1) then              
        if(numq>1) then              
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
        endif
     elseif(ioppq==1.and.ioprs==0) then             
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==0) then              
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
     endif
    elseif(icpqn==0.and.icrsn==1) then 
     if(ioppq==1.and.ioprs==1) then
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==1) then
        if(numq>1) then
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
        endif
     elseif(ioppq==1.and.ioprs==0) then
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==0) then
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
     endif
    endif

  !3333333333-----------------------------------------------------------
  elseif(ioppqn==1.and.ioprsn==0) then
    if(icpqn==0.and.icrsn==0) then 
     if(ioppq==1.and.ioprs==1) then
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==1) then              
        if(numq>1) then              
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
        endif
     elseif(ioppq==1.and.ioprs==0) then             
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==0) then              
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
            iin(idxn) = iout(idx)
        enddo
        enddo
        enddo
        enddo
     endif
   elseif(icpqn==1.and.icrsn==0) then
     if(ioppq==1.and.ioprs==1) then
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==1) then
        if(numq>1) then
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
        endif
     elseif(ioppq==1.and.ioprs==0) then
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
     elseif(ioppq==0.and.ioprs==0) then
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
            iin(idxn) = -iout(idx)
        enddo
        enddo
        enddo
        enddo
     endif
    endif
  !1111111111-----------------------------------------------------------
  elseif(ioppqn==0.and.ioprsn==0) then
         if(ioppq==1.and.ioprs==1) then
            do s = 1, nums
            do r = 1, numr
            do q = 1, numq
            do p = 1, nump
               idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                               + (q-1)*nump + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
               iin(idxn) = iout(idx)
            enddo
            enddo
            enddo
            enddo
         elseif(ioppq==0.and.ioprs==1) then              
            if(numq>1) then              
            do s = 1, nums
            do r = 1, numr
            do q = 2, numq
            do p = 1, q-1
               idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                               + (q-1)*(q-2)/2 + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
               iin(idxn) = iout(idx)
            enddo
            enddo
            enddo
            enddo
            endif
         elseif(ioppq==1.and.ioprs==0) then             
            do s = 2, nums
            do r = 1, s-1
            do q = 1, numq
            do p = 1, nump
               idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                               + (q-1)*nump + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
               iin(idxn) = iout(idx)
            enddo
            enddo
            enddo
            enddo
         elseif(ioppq==0.and.ioprs==0) then              
            do s = 2, nums
            do r = 1, s-1
            do q = 2, numq
            do p = 1, q-1
               idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                               + (q-1)*(q-2)/2 + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
!               write(6,*)'idx, idxn', idx, idxn
               iin(idxn) = iout(idx)
            enddo
            enddo
            enddo
            enddo
         endif
   endif
   return
   end

 
 


 
