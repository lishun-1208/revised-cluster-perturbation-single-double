! 2026.10.21
!> need initial dim, cyc, conv before use Davidson
module mod_davidson
    implicit none
    private

    abstract interface
        subroutine procedure_calAR(iCore, maxCore, istart, iAR, iR, symR)
            real*8, dimension(maxCore), intent(inout) :: iCore !> full work space
            integer*8, intent(in) :: maxCore !> size of iCore
            integer*8, intent(in) :: istart !> index of free space
            integer*8, intent(in) :: iAR !> index of AR
            integer*8, intent(in) :: iR !> index of R
            integer*8, intent(in) :: symR !> symmetry of R
        end subroutine procedure_calAR
    end interface

    public :: Davidson, procedure_calAR

contains
    ! write by lishun, 2024.9.25
    ! modified at 2026.4.29
    ! run Davidson iterations to solve AR=wR, dimA=n
    !-------------------------------------------------------------
    ! if want kth eiegen value of A:
    ! 1. give b_1,...,b_l as test vector, l>=k
    !    cal Ab_i, i=1:l
    ! 2. cal Abar_ij=(b_i,A*b_j)
    ! 3. solve Abar*r^(k)=a^(k)*r(k)
    ! 4. cal q^(k)=sum(i=1,l){r^(k)_i*b_i}
    ! 5. test convergence e=A*q^(k)-a^(k)*q^(k)
    ! 6. cal p_i=e_i/(a^(k)-A_ii), i=1:n, and orthonormalize to b
    !    then let b_(l+1) = p
    ! 7. return to step 1
    ! 8. if converged, let w=a^(k), R=q^(k)
    !---------------------------------------------------------------
    ! need subroutine: calAR(iCore, maxCore, istart, iAR, iR, symR)
    ! ATTENTATION: for non-Hermit A, the Abar is not Hermit, so eigenvector not ortongnal to each
    ! ATTENTATION: only support iop=0
    subroutine Davidson(iCore, maxCore, istart, fileEigenVV, &
                        diagA, dimA, nroot, istate, symR, calAR)
        use mod_common
        implicit none

        real*8, dimension(maxCore), intent(inout) :: iCore !> full work space
        integer*8, intent(in) :: maxCore !> size of iCore
        integer*8, intent(in) :: istart !> position of Free iCore
        ! file store initual guess and output eigenvalue and eigen vector for symR
        ! need with structure (val1,vec1),(val2,vec2),...
        character(len=*), intent(in) :: fileEigenVV
        real*8, dimension(dimA), intent(in) :: diagA !> diag element of A belong to a symmetry
        integer*8, intent(in) :: dimA !> length of eigenvector
        integer*8, intent(in) :: nroot !> total number of root need be calculate
        integer*8, intent(in) :: istate !> the first root need be calculate in all root of symR
        ! ie. if symR have total m root,
        !     then we calculate (istate,...,istate+nroot-1)th root
        integer*8, intent(in) :: symR !> symmetry of eiegnvector R
        procedure(procedure_calAR) :: calAR

        integer*8 :: iterations
        integer*8 :: dimAbar !> max dimension of Abar in step 2
        integer*8 :: dimIter !> dimension of Abar for an iteration
        integer*8 :: dimFree !> dimAbar-dimIter
        integer*8 :: ndimIter, ndimFree ! for new file when delet basis
        logical :: dim_exceed
        integer*8 :: dimDel !> the dimension will be delet when dim exceed, but need dimAbar-dimDel>nroot
        integer*8 :: iroot
        real*8, dimension(nroot) :: eigenVal !> store eigen value
        integer*8, dimension(nroot) :: rootConvCtrl !> 1 for convegenced, 0 for not
        real*8, dimension(:, :), allocatable :: Abar !> Abar in step 2
        real*8, dimension(:, :), allocatable :: tmpAbar
        ! eigenvectors and eigenvalues of Abar
        real*8, dimension(:), allocatable :: valBarR, valBarI ! eigenval of Abar
        real*8, dimension(:, :), allocatable :: vecBarR ! eigenvector of Abar
        real*8, dimension(:, :), allocatable :: vecNeed ! eigenvectors we need of Abar
        real*8, dimension(:, :), allocatable :: oldEigen, newEigen ! use when dim exceed, store eigenvector of Abar
        real*8, dimension(:), allocatable :: work ! work space for dgeev
        integer*8 :: lwork ! length of work of dgeev
        integer*8 :: info ! information of dgeev
        integer*8, dimension(:), allocatable :: lowestE
        real*8 :: tmp1, tmp2
        integer*8 :: nNewBase ! number new bases after an iteration
        integer*8 :: i, j, idx
        integer*8 :: i0, i1, i2, i3, nSize
        integer*8 :: maxSize
        real*8 :: maxError ! convergence contral
        character(len=10) :: file_b, file_Ab, nfile_b, nfile_Ab
        character(len=10) :: file_b1, file_Ab1, file_b2, file_Ab2
        real*8 :: tt0, tt1
        logical :: jump
        integer*8, external :: idsymsoc, irpdso, idamax
        real*8, external :: ddot, dnrm2

        if (nroot > 1) then
            write (6, '(A)') 'WARN: Davidson may have some error for nroot>1'
        end if

        if (iop /= 0) then
            write (6, '(A)') 'ERROR: Davidson only support iop=0'
            stop
        end if

        dimAbar = davidson_dim
        dimDel = dimAbar - nroot - 5
        if (nroot >= dimAbar) then
            write (6, *) 'ERROR: too much roots required', nroot, dimAbar
            return
        else if (nroot == 0) then
            return
        end if

        maxError = 10.d0**(-davidson_conv)
        write (6, *) 'Davidson convergence criteria:', maxError
        write (6, '(A,I)') 'max cycle for Davidson iteration:', davidson_cyc
        write (6, '(A,I)') 'max dimension of Davidson iteration matrix', davidson_dim

        eigenVal = 0.0d0
        ! max dimension of Abar, each iteration generate maximum nroot new vectors
        allocate (Abar(dimAbar, dimAbar))
        allocate (tmpAbar(dimAbar, dimAbar))
        allocate (valBarR(dimAbar))
        allocate (valBarI(dimAbar))
        allocate (vecBarR(dimAbar, dimAbar))
        allocate (vecNeed(dimAbar, nroot))
        vecNeed = 0.0d0
        allocate (lowestE(dimAbar))
        lwork = 10*dimAbar
        allocate (work(lwork))

        i0 = istart

        ! create file to store tmp data Ab and b in Davidson iteration
        file_b1 = 'davidson11'
        call IO_stream(0.0d0, iCore(i0), 1, 0, file_b1, '22', 'replace')
        file_Ab1 = 'davidson12'
        call IO_stream(0.0d0, iCore(i0), 1, 0, file_Ab1, '22', 'replace')
        file_b2 = 'davidson21'
        call IO_stream(0.0d0, iCore(i0), 1, 0, file_b2, '22', 'replace')
        file_Ab2 = 'davidson22'
        call IO_stream(0.0d0, iCore(i0), 1, 0, file_Ab2, '22', 'replace')
        ! file will use during iteration
        file_b = file_b1
        file_Ab = file_Ab1
        ! new file use when dimension exceed
        nfile_b = file_b2
        nfile_Ab = file_Ab2

        !------------- for initual guess
        call cpu_time(tt0)
        if (debug) then
            write (6, '(A)') 'start initial Davidson iteration'
        end if
        call initDavidson(iCore, maxCore, i0, dimA, dimAbar, nroot, istate, &
                          fileEigenVV, file_b, file_Ab, Abar, eigenVal, vecNeed, symR, &
                          calAR)
        call cpu_time(tt1)
        if (debug) then
            write (6, '(A,F)') 'CPU time for initial Davidson iteration', tt1 - tt0
        end if

        !--------- start iteration
        iterations = 1
        dimIter = nroot
        dimFree = dimAbar - dimIter
        dim_exceed = .false.
        iter: do
            call cpu_time(tt0)
            if (iterations > davidson_cyc) then
                write (6, *) 'WARN: Davidson iteration not convergence for all root'

                do iroot = 1, nroot
                if (rootConvCtrl(iroot) == 0) then
                    write (6, '(A,I3,A)') 'root', iroot, ' not covnverge'
                end if
                end do

                exit
            end if
            write (6, '(A,I5)') '--------- Davidson iteration:', iterations

            !-------------------------------------------------------
            ! solve Abar*r^(k)=a^(k)*r(k) in step 3
            !-------------------------------------------------------
            ! copy Abar, because dgeev will overwrite Abar
            tmpAbar = Abar
            call dgeev('N', 'V', dimIter, tmpAbar, dimAbar, valBarR, valBarI, &
                       vecBarR, dimAbar, vecBarR, dimAbar, work, lwork, info)
            if (info /= 0) then
                write (6, *) 'ERROR: dgeev return with error info:', info
                stop
            end if

            do i = 1, dimIter
                !---------- test whether eigenvalue be real, as we use only real part
                if (abs(valBarI(i)) > 1.0d-10) then
                    write (6, '(A,I4,2F)') 'WARN: eigenvalue be imag', i, valBarR(i), valBarI(i)
                end if

                !------ test eigensystem
                do j = 1, dimIter
                    ! A(j,k)*R(k,i)-w(i)*R(j,i)
                    iCore(i0 - 1 + j) = -valBarR(i)*vecBarR(j, i)
                    do idx = 1, dimIter
                        iCore(i0 - 1 + j) = iCore(i0 - 1 + j) + Abar(j, idx)*vecBarR(idx, i)
                    end do
                end do
                tmp1 = dnrm2(dimIter, iCore(i0), 1)
                if (abs(tmp1) > 1d-10) then
                    write (6, '(A,I4,F)') 'WARN: real part not eigenvector of Abar', i, tmp1
                end if

                !--------- test normalization
                tmp1 = dnrm2(dimIter, vecBarR(1, i), 1)
                if (abs(tmp1 - 1.0d0) > 1.0d-10) then
                    write (6, '(A,I5,F)') 'WARN: real part eigenvector of Abar not be normalized', &
                        i, tmp1
                    call dscal(dimIter, 1.0d0/tmp1, vecBarR(1, i), 1)
                end if
            end do

            !-------------------------------------------
            ! determinate the eigenvectors we need
            !-------------------------------------------
            !------------ sort eigenvalues
            call sort(dimIter, valBarR, lowestE)

            if (.not. lowest) then
                !------ determinate the eigenvectors closet to eigenvectors we need
                !------ just neet the (iroot)th element closet to 1
                ! I set eigenvalue I need to the lowest for convenient
                do iroot = 1, nroot
                    tmp1 = 0.0d0
                    ! cycle for vector
                    i = 1
                    do j = 1, dimIter
                        tmp2 = abs(vecBarR(iroot, lowestE(j))) ! vec_j(iroot)
                        if (tmp2 > tmp1) then
                            tmp1 = tmp2
                            i = j ! set vec_j as root(iroot)
                        end if
                    end do
                    idx = lowestE(i)

                    !------ resort lowestE
                    ! (iroot,...,i,...)->(i,irrot,...,...)
                    do j = 1, i - iroot
                        lowestE(i - j + 1) = lowestE(i - j)
                    end do
                    lowestE(iroot) = idx
                end do
            end if

            !------- determinate the eigenvectors use lowestE
            do iroot = 1, nroot
                ! copy vecBar to vecNeed
                vecNeed(:, iroot) = vecBarR(:, lowestE(iroot))
                eigenVal(iroot) = valBarR(lowestE(iroot))
                write (6, '(A,I4,A,F)') 'eigenvalue for state', iroot, ":", eigenVal(iroot)
            end do

            !-------------------------------------------------
            ! if dimension of Abar exceed up limit
            ! 1.delet the higest basis and store other basis in new file
            ! 2.generate new initial Abar
            ! in step 1, the eigenvector may not orthognal to each other,
            ! so we need find the transform to orthognal basis
            !-------------------------------------------------
            if (dim_exceed) then
                !--------- get orthonormalized eigenvector
                allocate (oldEigen(dimIter, ndimIter))
                allocate (newEigen(dimIter, ndimIter))
                if (debug) then
                    write (6, *) 'old dimIter=', dimIter
                end if
                ! get old vector
                do i = 1, ndimIter
                    !call dcopy(dimIter, vecBarR(1 + dimIter*(lowestE(i) - 1)), 1, &
                    !           oldEigen(1 + dimIter*(i - 1)), 1)
                    oldEigen(:, i) = vecBarR(:, lowestE(i))
                end do
                call Schmidt_2D(oldEigen, newEigen, dimIter, ndimIter)

                !--------- get new basis
                i1 = i0 + dimA
                i2 = i1 + dimA
                do iroot = 1, ndimIter
                    call dscal(dimA, 0.0d0, iCore(i0), 1)
                    call dscal(dimA, 0.0d0, iCore(i1), 1)
                    do i = 1, dimIter
                        tmp1 = newEigen(i, iroot) ! r^(k)_i
                        !----- calculate Aq, i0
                        ! read in Ab_i, i2
                        call IO_stream(0.0d0, iCore(i2), (i - 1)*dimA + 1, &
                                       dimA, file_Ab, '13', 'old')
                        ! Aq=r^(k)_i*Ab_i
                        call daxpy(dimA, tmp1, iCore(i2), 1, iCore(i0), 1)
                        !----- calculate q, i1
                        ! read in b_i, i2
                        call IO_stream(0.0d0, iCore(i2), (i - 1)*dimA + 1, &
                                       dimA, file_b, '13', 'old')
                        ! q = r^(k)_i*b_i
                        call daxpy(dimA, tmp1, iCore(i2), 1, iCore(i1), 1)
                    end do

                    !------ write Aq and q in new file
                    call IO_stream(0.0d0, iCore(i0), (iroot - 1)*dimA + 1, &
                                   dimA, nfile_Ab, '23', 'old')
                    call IO_stream(0.0d0, iCore(i1), (iroot - 1)*dimA + 1, &
                                   dimA, nfile_b, '23', 'old')
                end do
                deallocate (oldEigen, newEigen)
                write (6, *) 'finish generate new basis file'

                !--------- re-initial some variable
                file_Ab = nfile_Ab
                file_b = nfile_b
                if (file_Ab == file_Ab1) then
                    nfile_Ab = file_Ab2
                    nfile_b = file_b2
                else
                    nfile_Ab = file_Ab1
                    nfile_b = file_b1
                end if
                dimFree = ndimFree
                dimIter = ndimIter

                !--------- get new initial Abar
                ! calculate Abar_ji=(b_j,Ab_i)
                do i = 1, dimIter
                    ! get Ab_i, i0
                    call IO_stream(0.0d0, iCore(i0), (i - 1)*dimA + 1, &
                                   dimA, file_Ab, '13', 'old')
                    do j = 1, dimIter
                        ! get b_j, i1
                        call IO_stream(0.0d0, iCore(i1), (j - 1)*dimA + 1, &
                                       dimA, file_b, '13', 'old')
                        tmp1 = ddot(dimA, iCore(i0), 1, iCore(i1), 1)
                        Abar(j, i) = tmp1
                    end do
                end do
            end if

            nNewBase = 0
            rootConvCtrl = 0
            root: do iroot = 1, nroot
                i1 = i0 + dimA
                i2 = i1 + dimA
                i3 = i2 + dimA
                !----------------------------------------------------------
                !            calculate vector q and Aq (step 4)
                !----------------------------------------------------------
                if (dim_exceed) then
                    !------ read in Aq, i0
                    call IO_stream(0.0d0, iCore(i0), (iroot - 1)*dimA + 1, &
                                   dimA, file_Ab, '13', 'old')
                    !------ read in q, i1
                    call IO_stream(0.0d0, iCore(i1), (iroot - 1)*dimA + 1, &
                                   dimA, file_b, '13', 'old')
                else
                    call dscal(dimA, 0.0d0, iCore(i0), 1)
                    call dscal(dimA, 0.0d0, iCore(i1), 1)
                    do i = 1, dimIter
                        tmp1 = vecNeed(i, iroot) ! r^(k)_i
                        !----- calculate Aq, i0
                        ! read in Ab_i, i2
                        call IO_stream(0.0d0, iCore(i2), (i - 1)*dimA + 1, &
                                       dimA, file_Ab, '13', 'old')
                        ! Aq=r^(k)_i*Ab_i
                        call daxpy(dimA, tmp1, iCore(i2), 1, iCore(i0), 1)
                        !----- calculate q, i1
                        ! read in b_i, i2
                        call IO_stream(0.0d0, iCore(i2), (i - 1)*dimA + 1, &
                                       dimA, file_b, '13', 'old')
                        ! q = r^(k)_i*b_i
                        call daxpy(dimA, tmp1, iCore(i2), 1, iCore(i1), 1)
                    end do
                end if

                tmp1 = dnrm2(dimA, iCore(i1), 1)
                if (abs(tmp1 - 1.0d0) > 1.0d-10) then
                    write (6, '(A,I4,F)') 'WARN: eigenvector of A not normzlized', &
                        iroot, tmp1
                    call dscal(dimA, 1/tmp1, iCore(i1), 1)
                    call dscal(dimA, 1/tmp1, iCore(i0), 1)
                end if

                !--------- write eigenvalue and eigenvector in hard disk
                idx = (istate + iroot - 2)*(1 + dimA) + 1
                call IO_stream(eigenVal(iroot), iCore(i1), idx, &
                               dimA, fileEigenVV, '21', 'old')

                !------------------------------------------------------------------------
                !       calculate error vector and test convergence (step 5)
                !------------------------------------------------------------------------
                !---------- calculate error vector e=Aq-aq, i0
                call daxpy(dimA, -eigenVal(iroot), iCore(i1), 1, iCore(i0), 1)

                !----------- test convergence
                tmp1 = iCore(i0 + idamax(dimA, iCore(i0), 1) - 1)
                write (6, '(A,I3,A,I4,A,F)') 'max element of error vector for root', &
                    iroot, ' after', iterations, ' iterations:', tmp1
                if (abs(tmp1) < maxError) then
                    write (6, '(A,I3,A)') 'Davidson for root', iroot, ' converged'
                    rootConvCtrl(iroot) = 1
                    cycle ! not calculate new base vector for converged root
                end if

                !------------------------------------------------------
                ! if not convergence, calculate new base vector (step 6)
                !------------------------------------------------------
                ! check dimension of Abar, not generate new base if exceed
                if (dimFree == 0) then
                    cycle
                end if

                !--------- calculate vector p_i=e_i/(a^(k)-A_ii), i1
                do i = 1, dimA
                    tmp1 = eigenVal(iroot) - diagA(i)
                    if (abs(tmp1) < 1.0d-4) then
                        write (6, *) 'WARN: divisor is too small', tmp1
                        if (tmp1 < 0.0d0) then
                            tmp1 = -1.0d-4 ! avoid divisor be 0
                        else
                            tmp1 = 1.0d-4 ! avoid divisor be 0
                        end if
                    end if
                    iCore(i1 - 1 + i) = iCore(i0 - 1 + i)/tmp1
                end do

                !--------- orthonormalize p
                jump = .false.
                call SchmidtS(iCore(i1), dimA, dimIter + nNewBase, &
                              file_b, iCore(i2), maxCore - i2 + 1, jump)
                if (jump) then
                    cycle
                end if
                !--------- write in new basis vector
                call IO_stream(0.0d0, iCore(i1), (dimIter + nNewBase)*dimA + 1, &
                               dimA, file_b, '23', 'old')

                !--------- calculate and write in Ap, i2
                call calAR(iCore, maxCore, i3, i2, i1, symR)
                call IO_stream(0.0d0, iCore(i2), (dimIter + nNewBase)*dimA + 1, &
                               dimA, file_Ab, '23', 'old')

                dimFree = dimFree - 1
                nNewBase = nNewBase + 1
            end do root

            !--------- all root converged if i=nroot
            i = 0
            do iroot = 1, nroot
                i = i + rootConvCtrl(iroot)
            end do
            if (i == nroot) then
                write (6, '(A)') 'All roots converged'
                exit
            end if

            !--------- if dimFree=0, next iteration will delet some base vector
            if (dim_exceed) then
                dim_exceed = .false.
            end if
            if (dimFree == 0) then
                dim_exceed = .true.
                if (debug) then
                    write (6, *) 'Dimension of Abar exceed up limit, next iteration will delet some basis'
                end if
                ndimIter = dimAbar - dimDel
                if (ndimIter < nroot) then
                if (debug) then
                    write (6, '(A)') 'ndimIter<nroot, set it be nroot'
                end if
                ndimIter = nroot
                ndimFree = dimAbar - ndimIter
                else
                ndimFree = dimDel
                end if
            end if

            !----------------------------------------------------------
            !                  calculate Abar (step 2)
            !----------------------------------------------------------
            !--------- calculate new elements Abar_ij = (b_i,Ab_j)
            i1 = i0 + dimA
            do j = 1, dimIter
                ! read in Ab_j, i0
                call IO_stream(0.0d0, iCore(i0), (j - 1)*dimA + 1, &
                               dimA, file_Ab, '13', 'old')
                do i = dimIter + 1, dimIter + nNewBase
                    ! read in b_i, i1
                    call IO_stream(0.0d0, iCore(i1), (i - 1)*dimA + 1, &
                                   dimA, file_b, '13', 'old')
                    tmp1 = ddot(dimA, iCore(i1), 1, iCore(i0), 1)
                    Abar(i, j) = tmp1
                end do
            end do
            do j = dimIter + 1, dimIter + nNewBase
                ! read in Ab_j, i0
                call IO_stream(0.0d0, iCore(i0), (j - 1)*dimA + 1, &
                               dimA, file_Ab, '13', 'old')
                do i = 1, dimIter + nNewBase
                    ! read in b_i, i1
                    call IO_stream(0.0d0, iCore(i1), (i - 1)*dimA + 1, &
                                   dimA, file_b, '13', 'old')
                    tmp1 = ddot(dimA, iCore(i1), 1, iCore(i0), 1)
                    Abar(i, j) = tmp1
                end do
            end do

            dimIter = dimIter + nNewBase
            iterations = iterations + 1
            call cpu_time(tt1)
            if (debug) then
                write (6, '(A,I4,A,F)') 'CPU time for iter', iterations, ':', tt1 - tt0
            end if
        end do iter

        deallocate (Abar, valBarR, valBarI, vecBarR, &
                    vecNeed, lowestE, work)

        !-----------------------------------
        ! test eigenvalue and eigenvector
        !-----------------------------------
        i1 = i0 + dimA
        i2 = i1 + dimA
        do iroot = 1, nroot
            !---- get R and w, i0
            idx = (istate + iroot - 2)*(1 + dimA) + 1
            call IO_stream(eigenVal(iroot), iCore(i0), idx, dimA, fileEigenVV, '11', 'old')
            !---- e = AR-wR, i1
            call calAR(iCore, maxCore, i2, i1, i0, symR)
            call daxpy(dimA, -eigenVal(iroot), iCore(i0), 1, iCore(i1), 1)

            i = idamax(dimA, iCore(i1), 1)
            tmp1 = iCore(i1 - 1 + i)
            if (abs(tmp1) > maxError) then
                write (6, '(A,I3,A,F)') 'WARN: root', iroot, ' not pass eigen test', tmp1
            else
                write (6, '(A,I3,A)') 'root', iroot, ' pass test'
            end if
        end do

        !------ test orthogonalize
        !i1 = i0 + dimA
        !i2 = i1 + dimA
        !do i = 1, dimIter
        !    call IO_stream(0.0d0, iCore(i0), dimA*(i - 1) + 1, dimA, file_b, '13', 'old')
        !    tmp1 = dnrm2(dimA, iCore(i0), 1)
        !    if (abs(tmp1 - 1.0d0) > 1.0d-10) then
        !        write (6, '(A,I3,F)') 'WARN: norm:', i, tmp1
        !    end if
        !    do j = i + 1, dimIter
        !        call IO_stream(0.0d0, iCore(i1), dimA*(j - 1) + 1, dimA, file_b, '13', 'old')
        !        tmp1 = ddot(dimA, iCore(i0), 1, iCore(i1), 1)
        !        if (abs(tmp1) > 1.0d-10) then
        !            write (6, '(A,2I3,F)') 'WARN: orth:', i, j, tmp1
        !        end if
        !    end do
        !end do
    end subroutine Davidson

    subroutine initDavidson(iCore, maxCore, istart, dimA, dimAbar, nroot, istate, &
                            fileEigenVV, file_b, file_Ab, Abar, eigenVal, vecNeed, symR, &
                            calAR)
        implicit none

        real*8, dimension(maxCore), intent(inout) :: iCore
        integer*8, intent(in) :: maxCore, istart
        integer*8, intent(in) :: dimA, dimAbar, nroot, istate
        character(len=*), intent(in) :: fileEigenVV
        character(len=10), intent(in) :: file_b, file_Ab
        real*8, dimension(dimAbar, dimAbar), intent(out) :: Abar
        real*8, dimension(nroot), intent(out) :: eigenVal
        real*8, dimension(dimAbar, nroot), intent(out) :: vecNeed
        integer*8, intent(in) :: symR
        procedure(procedure_calAR) :: calAR

        integer*8 :: i0, i1, i2, i, j, idx
        real*8, external :: dnrm2, ddot

        i0 = istart
        i1 = i0 + dimA
        i2 = i1 + dimA
        !--------- cycle for root
        do i = 1, nroot
            ! read in initial guess b_i, i0
            call IO_stream(eigenVal(i), iCore(i0), (istate + i - 2)*(dimA + 1) + 1, &
                           dimA, fileEigenVV, '11', 'old')
            ! write in b_i
            call IO_stream(0.0d0, iCore(i0), (i - 1)*dimA + 1, dimA, file_b, '23', 'old')

            ! calculate Ab_i, i1
            call calAR(iCore, maxCore, i2, i1, i0, symR)
            ! write in Ab_i
            call IO_stream(0.0d0, iCore(i1), (i - 1)*dimA + 1, dimA, file_Ab, '23', 'old')

            ! calculate Abar_ji=(b_j,Ab_i)
            do j = 1, nroot
                ! read in b_j, i0
                call IO_stream(0.0d0, iCore(i0), (istate + j - 2)*(dimA + 1) + 2, &
                               dimA, fileEigenVV, '13', 'old')
                Abar(j, i) = ddot(dimA, iCore(i0), 1, iCore(i1), 1)
            end do

            ! get initial vecNeed
            vecNeed(i, i) = 1.0d0
        end do

    end subroutine initDavidson

    ! write by lishun, 2024.9.16
    ! sort element in vector from lower to higher
    subroutine sort(n, vecin, idx)
        implicit none

        integer*8, intent(in) :: n ! dimension of vecin and idx
        real*8, dimension(n), intent(in) :: vecin ! vector need be sort
        integer*8, dimension(n), intent(out) :: idx ! idx from lower to higher

        integer*8 :: i, j
        integer*8 :: tmp

        idx = (/(i, i=1, n)/)
        do i = 1, n - 1
            do j = i + 1, n
                if (vecin(idx(j)) < vecin(idx(i))) then
                    tmp = idx(j)
                    idx(j) = idx(i)
                    idx(i) = tmp
                end if
            end do
        end do

    end subroutine sort

    ! write at 2025.5.27
    ! 2D Schmidt
    ! need the old basis be normalized
    subroutine Schmidt_2D(inM, outM, length, dim)
        implicit none

        real*8, dimension(length, dim), intent(in) :: inM
        real*8, dimension(length, dim), intent(out) :: outM
        integer*8, intent(in) :: length, dim

        integer*8 :: i, j
        real*8 :: norm, overlap
        real*8, external :: dnrm2, ddot

        outM = inM

        !--------- for other vectors
        do i = 2, dim ! new basis
            !------ orthognalize to new basis
            do j = 1, i - 1 ! new basis
                ! |ni>=|ni>-<nj|ni>*|nj>
                overlap = ddot(length, outM(1, j), 1, outM(1, i), 1)
                call daxpy(length, -overlap, outM(1, j), 1, outM(1, i), 1)
            end do
            !------ normalize |i>
            norm = dnrm2(length, outM(1, i), 1)
            call dscal(length, 1/norm, outM(1, i), 1)
        end do

        !--------- test
        do i = 1, dim
            norm = dnrm2(length, outM(1, i), 1)
            if (abs(norm - 1.0d0) > 1.0d-14) then
                write (6, *) 'WARN: Schmidt not normalized:', i, norm
            end if
            do j = i + 1, dim
                overlap = ddot(length, outM(1, j), 1, outM(1, i), 1)
                if (abs(overlap) > 1.0d-14) then
                    write (6, *) 'WARN: Schmidt not orthogonalized:', i, j, overlap
                end if
            end do
        end do
    end subroutine Schmidt_2D

    ! write by lishun, 2024.9.27
    ! run Schmidt orthonormalize procedure for vecin
    ! other vectors must be orthonormalized, ie. <i|i>=1
    ! only single
    !----------------------------------
    ! modified at 2025.5.27
    subroutine SchmidtS(vecin, length, dim, file_b, iCore, maxCore, jump)
        implicit none

        real*8, dimension(length), intent(inout) :: vecin ! |vec>
        integer*8, intent(in) :: length ! length of vector
        integer*8, intent(in) :: dim ! number of vectors which vecin need be orthogonal
        character(len=10), intent(in) :: file_b
        real*8, dimension(maxCore), intent(inout) :: iCore
        integer*8, intent(in) :: maxCore
        logical, intent(inout) :: jump

        integer*8 :: i, j
        integer*8 :: i0
        real*8 :: overlop, norm
        real*8, external :: ddot, dnrm2

        i0 = 1
        ! normalize |vec>
        norm = dnrm2(length, vecin, 1)
        call dscal(length, 1/norm, vecin, 1)

        do i = 1, dim
            ! read in |i>, i0
            call IO_stream(0.0d0, iCore(i0), (i - 1)*length + 1, length, file_b, '13', 'old')
            ! calculate overlop=<i|vec>
            overlop = ddot(length, iCore(i0), 1, vecin, 1)
            ! orthogonal to |i>, |vec>=|vec>-<i|vec>*|i>
            call daxpy(length, -overlop, iCore(i0), 1, vecin, 1)
        end do

        ! normalize |vec>
        norm = dnrm2(length, vecin, 1)
        if (norm < 1.0d-6) then
            write (6, *) 'WARN: new vector is linear independent:', norm
            jump = .true.
            return
        end if
        call dscal(length, 1/norm, vecin, 1)

    end subroutine SchmidtS

end module mod_davidson
