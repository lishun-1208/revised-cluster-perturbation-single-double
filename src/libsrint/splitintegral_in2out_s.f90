  subroutine splitintegral_in2out_s(nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, ioffspin, iout, &
                   numpn, numqn, ioppqn, icpqn, numrn, numsn, ioprsn, icrsn, numpqn, numrsn, & 
                   ioffp, ioffq, ioffr, ioffs, ioffnspin, iin)
  implicit none 
  integer,intent(in):: nump, numq, ioppq, numr, nums, ioprs, numpq, numrs, ioffspin
  integer,intent(in):: numpn, numqn, ioppqn, numrn, numsn, ioprsn, numpqn, numrsn 
  integer, intent(in):: ioffp, ioffq, ioffr, ioffs, ioffnspin, icpqn, icrsn
  real*4, intent(in) :: iin(*)
  real*4, intent(out) :: iout(*)
  integer p, q, r, s, idx, idxn
 ! write(6,*)'iop', ioppqn, ioprsn, ioppq, ioprs
  !1111111111-----------------------------------------------------------
  if(ioppqn==1.and.ioprsn==1) then 
  !1222222222-----------------------------------------------------------
    if(icpqn==0.and.icrsn==0) then 
      if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==1) then 
         if(numq>1) then              
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
         endif
      elseif(ioppq==1.and.ioprs==0) then             
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==0) then              
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      endif
  !1122222222-----------------------------------------------------------
    elseif(icpqn==0.and.icrsn==1) then 
      if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==1) then
         if(numq>1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
         endif
      elseif(ioppq==1.and.ioprs==0) then
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==0) then
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffq+q-1)*numpn + ioffp+p
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      endif
  !1133333333-----------------------------------------------------------
    elseif(icpqn==1.and.icrsn==0) then 
      if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==1) then
         if(numq>1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
         endif
      elseif(ioppq==1.and.ioprs==0) then
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==0) then
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = -iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      endif
  !1144444444-----------------------------------------------------------
    elseif(icpqn==1.and.icrsn==1) then 
      if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*numr+r-1)*numpq &
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==1) then
         if(numq>1) then
!$OMP PARALLEL DO
         do s = 1, nums
         do r = 1, numr
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
         endif
      elseif(ioppq==1.and.ioprs==0) then
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 1, numq
         do p = 1, nump
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*nump + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      elseif(ioppq==0.and.ioprs==0) then
!$OMP PARALLEL DO
         do s = 2, nums
         do r = 1, s-1
         do q = 2, numq
         do p = 1, q-1
            idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                            + (q-1)*(q-2)/2 + p
            idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                              + (ioffp+p-1)*numqn + ioffq+q
            iout(idx) = iin(idxn)
         enddo
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      endif
     endif

  !2222222222-----------------------------------------------------------
  elseif(ioppqn==0.and.ioprsn==1) then
    if(icpqn==0.and.icrsn==0) then  
     if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==1) then              
        if(numq>1) then              
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
        endif
     elseif(ioppq==1.and.ioprs==0) then             
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==0) then              
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*numrn+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     endif
    elseif(icpqn==0.and.icrsn==1) then 
     if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==1) then
        if(numq>1) then
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
        endif
     elseif(ioppq==1.and.ioprs==0) then
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==0) then
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffr+r-1)*numsn+ioffs+s-1)*numpqn&
                             + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     endif
    endif

  !3333333333-----------------------------------------------------------
  elseif(ioppqn==1.and.ioprsn==0) then
    if(icpqn==0.and.icrsn==0) then 
     if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==1) then              
        if(numq>1) then              
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
        endif
     elseif(ioppq==1.and.ioprs==0) then             
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==0) then              
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffq+q-1)*numpn + ioffp+p
           iout(idx) = iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     endif
   elseif(icpqn==1.and.icrsn==0) then
     if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==1) then
        if(numq>1) then
!$OMP PARALLEL DO
        do s = 1, nums
        do r = 1, numr
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
        endif
     elseif(ioppq==1.and.ioprs==0) then
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 1, numq
        do p = 1, nump
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*nump + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     elseif(ioppq==0.and.ioprs==0) then
!$OMP PARALLEL DO
        do s = 2, nums
        do r = 1, s-1
        do q = 2, numq
        do p = 1, q-1
           idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                           + (q-1)*(q-2)/2 + p
           idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                             + (ioffp+p-1)*numqn + ioffq+q
           iout(idx) = -iin(idxn)
        enddo
        enddo
        enddo
        enddo
!$OMP END PARALLEL DO
     endif
    endif
  !1111111111-----------------------------------------------------------
  elseif(ioppqn==0.and.ioprsn==0) then
         if(ioppq==1.and.ioprs==1) then
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, numr
            do q = 1, numq
            do p = 1, nump
               idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                               + (q-1)*nump + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
               iout(idx) = iin(idxn)
            enddo
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
         elseif(ioppq==0.and.ioprs==1) then              
            if(numq>1) then              
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, numr
            do q = 2, numq
            do p = 1, q-1
               idx = ioffspin + ((s-1)*numr+r-1)*numpq&
                               + (q-1)*(q-2)/2 + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
               iout(idx) = iin(idxn)
            enddo
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
            endif
         elseif(ioppq==1.and.ioprs==0) then             
!$OMP PARALLEL DO
            do s = 2, nums
            do r = 1, s-1
            do q = 1, numq
            do p = 1, nump
               idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                               + (q-1)*nump + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
               iout(idx) = iin(idxn)
            enddo
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
         elseif(ioppq==0.and.ioprs==0) then              
!$OMP PARALLEL DO
            do s = 2, nums
            do r = 1, s-1
            do q = 2, numq
            do p = 1, q-1
               idx = ioffspin + ((s-1)*(s-2)/2+r-1)*numpq&
                               + (q-1)*(q-2)/2 + p
               idxn = ioffnspin + ((ioffs+s-1)*(ioffs+s-2)/2+ioffr+r-1)*numpqn&
                                 + (ioffq+q-1)*(ioffq+q-2)/2 + ioffp+p
!               write(6,*)'idx, idxn', idx, idxn
                         iout(idx) = iin(idxn)
            enddo
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
         endif
   endif
   return
   end

 
 


 
