! lish, 2026.7.8
module mod_operator
    implicit none
    private

    logical, parameter, public :: T = .true., F = .false.

    type, public :: opt
        ! operator I(pq,rs)
        integer*8, dimension(8) :: popp, popq, popr, pops ! population of p, q, r, and s belong to irrep
        integer*8 :: ioppq, ioprs ! ioppq=0 for p<q; 1 for p, q; 2 for p<=q
        integer*8 :: sym ! symmetry of operator I
        integer*8 :: idx ! the position of operator I in work space iCore
        integer*8 :: len ! length of I(pq,rs)
        character(len=8) :: f ! the file which store I
        logical :: incore ! whether store in memory
        real*8, dimension(:), pointer :: arr ! the pointer to the array
    contains
        ! key word pass will make operator be firse parameter in procerdure
        procedure, pass :: init => init_opt
        procedure, pass :: print => print_opt
        procedure, pass :: norm => norm_opt
        procedure, pass :: arr_copy => arr_copy_opt
        procedure, pass :: arr_clear => arr_clear_opt !> set 0 for arr
    end type opt

    type, public :: SAopt
        integer*8, dimension(8) :: poppa, poppb
        integer*8, dimension(8) :: popqa, popqb
        integer*8 :: ioppq
        integer*8, dimension(8) :: popra, poprb
        integer*8, dimension(8) :: popsa, popsb
        integer*8 :: ioprs
        integer*8 :: sym
        integer*8 :: particle !> number of particles
        logical :: incore
        integer*8 :: size !> total size in memory
        logical :: UHF !> reference state
        type(opt) :: AA, BB, AB, BA
        logical, dimension(4) :: single_exist = F !> exist of AA, BB, AB, BA part
        type(opt) :: AAAA, AABB, AAAB, AABA
        type(opt) :: BBAA, BBBB, BBAB, BBBA
        type(opt) :: ABAA, ABBB, ABAB, ABBA
        type(opt) :: BAAA, BABB, BAAB, BABA
        logical, dimension(4, 4) :: double_exist = F !> exist of AA, BB, AB, BA part
        logical :: SOC = F
    contains
        procedure, pass :: init1_SAopt !> initialize 1 particle
        procedure, pass :: init2_SAopt !> initialize 2 particle
        generic :: init => init1_SAopt, init2_SAopt
        procedure, pass :: norm => norm_SAopt !> print norm for each spin part
        procedure, pass :: print => print_SAopt !> print information
        !procedure, pass :: spin_exist_SAopt1 !> test whether require spin 1 exist
        !procedure, pass :: spin_exist_SAopt2 !> test whether require spin 2 exist
        !generic :: spin_exist => spin_exist_SAopt1, spin_exist_SAopt2
        procedure, pass :: spin_exist => spin_exist_SAopt !> test whether require spin exist
        procedure, pass :: arr_copy => arr_copy_SAopt !> copy array for all exist spin
        procedure, pass :: arr_clear => arr_clear_SAopt !> set 0 for all exist spin
        procedure, pass :: copy => copy_SAopt !> copy object but not include array and file name, only incore
        procedure, pass :: read => read_SAopt !> getallso for all spin case
        procedure, pass :: write => write_SAopt !> putallso for all spin case
        procedure, pass :: create_f => create_f_SAopt !> clmoioso0 for all spin case
        procedure, pass :: add => add_SAopt !> self += fact*in
        procedure, pass :: transpose => transpose_SAopt !> itranspso for all spin case
    end type SAopt

    public :: copy_opt, pop_eq

contains
    subroutine init_opt(V, iCore, popp, popq, popr, pops, ioppq, ioprs, sym, &
                        incore, idx, file)
        implicit none

        class(opt), intent(inout) :: V
        real*8, dimension(*), intent(in) :: iCore !> the full work space
        integer*8, dimension(8), intent(in) :: popp, popq, popr, pops ! population of p, q, r, and s belong to irrep
        integer*8, intent(in) :: ioppq, ioprs ! ioppq=0 for p<q; 1 for p, q; 2 for p<=q
        integer*8, intent(in) :: sym !> symmetry of operator I
        logical, intent(in), optional :: incore !> whether store in memory, default=T
        integer*8, intent(in), optional :: idx !> the position of operator I in iCore
        character(len=8), intent(in), optional :: file !> the file which store I

        target :: iCore
        integer*8, external :: idsymsoc

        V%popp = popp
        V%popq = popq
        V%popr = popr
        V%pops = pops
        V%ioppq = ioppq
        V%ioprs = ioprs
        V%sym = sym
        V%len = idsymsoc(V%sym, V%popp, V%popq, V%ioppq, V%popr, V%pops, V%ioprs)

        if (.not. present(incore)) then
            V%incore = T
        else
            V%incore = incore
        end if

        if (V%incore) then
            V%idx = idx
            V%arr => iCore(V%idx:V%idx)
        end if

        if ((.not. V%incore) .and. (.not. present(file))) then
            write (6, *) 'ERROR: init opt class error'
            stop
        end if

        if (present(file)) then
            V%f = file
        end if

    end subroutine init_opt

    subroutine print_opt(V)
        implicit none

        class(opt), intent(in) :: V

        integer*8 :: i
        real*8 :: norm

        write (6, '(A,8I5)') 'popp: ', (V%popp(i), i=1, 8)
        write (6, '(A,8I5)') 'popq: ', (V%popq(i), i=1, 8)
        write (6, '(A,8I5)') 'popr: ', (V%popr(i), i=1, 8)
        write (6, '(A,8I5)') 'pops: ', (V%pops(i), i=1, 8)
        write (6, '(A,I1)') 'ioppq: ', V%ioppq
        write (6, '(A,I1)') 'ioprs: ', V%ioprs
        write (6, '(A,I1)') 'symmetry: ', V%sym
        write (6, '(A,I)') 'index: ', V%idx
        write (6, '(A,I)') 'length: ', V%len
        write (6, '(2A)') 'file: ', V%f
        write (6, '(A,L)') 'store in memory: ', V%incore

    end subroutine print_opt

    !> calculate norm of an I(pq,rs)
    !> ATTENTION: the parentheses is necessary
    function norm_opt(V, iCore, maxCore) result(norm)
        implicit none

        class(opt), intent(in) :: V
        real*8, dimension(*), intent(inout), optional :: iCore !> free work space
        integer*8, intent(in), optional :: maxCore
        real*8 :: norm

        real*8, external :: dnrm2

        if (V%incore) then
            norm = dnrm2(V%len, V%arr, 1)
        else
            call calLargeNorm(norm, V%f, V%popp, V%popq, V%ioppq, &
                              V%popr, V%pops, V%ioprs, V%sym, iCore, maxCore)
        end if

    end function norm_opt

    !> copy an array from another array
    subroutine arr_copy_opt(opt_out, opt_in)
        implicit none

        class(opt), intent(inout) :: opt_out
        type(opt), intent(in) :: opt_in

        if (opt_in%len /= opt_out%len) then
            write (6, *) 'ERROR: arr_copy_opt len not be the same'
            stop
        end if

        call dcopy(opt_in%len, opt_in%arr, 1, opt_out%arr, 1)

    end subroutine arr_copy_opt

    subroutine arr_clear_opt(self)
        implicit none

        class(opt), intent(inout) :: self

        call dscal(self%len, 0.0d0, self%arr, 1)

    end subroutine arr_clear_opt

    !> copy a opt obejct to an other, except the position
    subroutine copy_opt(in, out, incore, iCore, idx, file)
        implicit none

        type(opt), intent(in) :: in
        type(opt), intent(out) :: out
        logical, intent(in), optional :: incore !> whether store in memory
        real*8, dimension(*), intent(in), optional :: iCore !> full work space
        integer*8, intent(in), optional :: idx !> position of arr
        character(len=8), intent(in), optional :: file

        target :: iCore

        out%popp = in%popp
        out%popq = in%popq
        out%popr = in%popr
        out%pops = in%pops
        out%ioppq = in%ioppq
        out%ioprs = in%ioprs
        out%sym = in%sym
        out%len = in%len

        if (.not. present(incore)) then
            out%incore = T
        else
            out%incore = incore
        end if

        if (out%incore) then
            out%idx = idx
            out%arr => iCore(out%idx:out%idx)
        end if

        if ((.not. out%incore) .and. (.not. present(file))) then
            write (6, *) 'ERROR: copy_opt error'
            stop
        end if

        if (present(file)) then
            out%f = file
        end if

    end subroutine copy_opt

    subroutine init1_SAopt(self, poppa, poppb, popqa, popqb, sym, ref, &
                           incore, iCore, idx, n_spin, spin, file)
        use mod_orbit
        implicit none

        ! sym and incore are global initialized
        class(SAopt), intent(out) :: self
        integer*8, dimension(8), intent(in) :: poppa, poppb
        integer*8, dimension(8), intent(in) :: popqa, popqb
        integer*8, intent(in) :: sym
        logical*8, intent(in) :: ref !> whether be UHF reference
        logical, intent(in), optional :: incore !> whether store in memory, default=T
        real*8, dimension(*), intent(in), optional :: iCore !> work space
        integer*8, intent(in), optional :: idx !> start position in memery
        ! if incore, will store as the following order
        integer*8, intent(in), optional :: n_spin !> number of spin case
        character(len=2), dimension(n_spin), intent(in), optional :: spin !> spin case
        character(len=8), dimension(n_spin), intent(in), optional :: file !> file

        integer*8 :: i, i0, i1
        character(len=4), dimension(4) :: spin_case
        integer*8 :: ioppq = 1
        character(len=8), dimension(4) :: spin_f

        self%poppa = poppa
        self%poppb = poppb
        self%popqa = popqa
        self%popqb = popqb
        self%popra = zorb
        self%poprb = zorb
        self%popsa = zorb
        self%popsb = zorb
        self%ioppq = ioppq
        self%ioprs = 1
        self%sym = sym
        self%particle = 1
        self%size = 0
        self%UHF = ref
        if (.not. present(incore)) then
            self%incore = T
        else
            self%incore = incore
        end if

        if (present(spin)) then
            spin_case(1:n_spin) = spin(1:n_spin)
        else
            spin_case = ['AA', 'BB', 'AB', 'BA']
        end if

        if (present(file)) then
            spin_f(1:n_spin) = file(1:n_spin)
        else
            if (.not. self%incore) then
                write (6, '(A)') 'ERROR: init1_SAopt outcore and not give file'
                stop
            end if
            spin_f = 'SPINFILE'
        end if

        i0 = idx
        ! cycle for spin case
        do i = 1, n_spin
            if (spin_case(i) == 'AA') then
                call self%AA%init(iCore, poppa, popqa, zorb, zorb, &
                                  ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%single_exist(1) = T
                if (self%incore) then
                    i0 = i0 + self%AA%len
                end if
            else if (spin_case(i) == 'BB') then
                call self%BB%init(iCore, poppb, popqb, zorb, zorb, &
                                  ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%single_exist(2) = T
                if (self%incore) then
                    i0 = i0 + self%BB%len
                end if
            else if (spin_case(i) == 'AB') then
                call self%AB%init(iCore, poppa, popqb, zorb, zorb, &
                                  ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%single_exist(3) = T
                if (self%incore) then
                    i0 = i0 + self%AB%len
                end if
            else if (spin_case(i) == 'BA') then
                call self%BA%init(iCore, poppb, popqa, zorb, zorb, &
                                  ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%single_exist(4) = T
                if (self%incore) then
                    i0 = i0 + self%BA%len
                end if
            else
                write (6, '(2A)') 'ERROR: spin case in init1_SAopt error.', spin_case(i)
                stop
            end if
        end do
        self%size = i0 - idx

    end subroutine init1_SAopt

    subroutine init2_SAopt(self, poppa, poppb, popqa, popqb, &
                           popra, poprb, popsa, popsb, ioppq, ioprs, sym, ref, &
                           incore, iCore, idx, n_spin, spin, file)
        implicit none

        ! iop, sym, and incore are global initialized
        class(SAopt), intent(out) :: self
        integer*8, dimension(8), intent(in) :: poppa, poppb
        integer*8, dimension(8), intent(in) :: popqa, popqb
        integer*8, dimension(8), intent(in) :: popra, poprb
        integer*8, dimension(8), intent(in) :: popsa, popsb
        integer*8, intent(in) :: ioppq, ioprs
        integer*8, intent(in) :: sym
        logical, intent(in) :: ref !> whether be UHF reference
        logical, intent(in), optional :: incore !> whether store in memory, default=T
        real*8, dimension(*), intent(in), optional :: iCore !> work space
        integer*8, intent(in), optional :: idx !> start position in memery
        ! if incore, will store as the following order
        integer*8, intent(in), optional :: n_spin !> number of spin case
        character(len=4), dimension(n_spin), intent(in), optional :: spin !> spin case
        character(len=8), dimension(n_spin), intent(in), optional :: file !> file

        integer*8 :: i, i0, i1
        character(len=4), dimension(16) :: spin_case
        logical :: inout
        character(len=8), dimension(16) :: spin_f

        self%poppa = poppa
        self%poppb = poppb
        self%popqa = popqa
        self%popqb = popqb
        self%popra = popra
        self%poprb = poprb
        self%popsa = popsa
        self%popsb = popsb
        self%ioppq = ioppq
        self%ioprs = ioprs
        self%sym = sym
        self%particle = 2
        self%size = 0
        self%UHF = ref
        if (.not. present(incore)) then
            self%incore = T
        else
            self%incore = incore
        end if

        if (ioppq /= 1) then
            if ((.not. pop_eq(poppa, popqa)) .or. (.not. pop_eq(poppb, popqb))) then
                write (6, *) 'ERRPR: init2_SAopt error with ioppq.'
                stop
            end if
        end if
        if (ioprs /= 1) then
            if ((.not. pop_eq(popra, popsa)) .or. (.not. pop_eq(poprb, popsb))) then
                write (6, *) 'ERRPR: init2_SAopt error with ioprs.'
                stop
            end if
        end if

        if (present(spin)) then
            spin_case(1:n_spin) = spin(1:n_spin)
        else
            spin_case = ['AAAA', 'AABB', 'AAAB', 'AABA', &
                         'BBAA', 'BBBB', 'BBAB', 'BBBA', &
                         'ABAA', 'ABBB', 'ABAB', 'ABBA', &
                         'BAAA', 'BABB', 'BAAB', 'BABA']
        end if

        if (present(file)) then
            spin_f(1:n_spin) = file(1:n_spin)
        else
            if (.not. self%incore) then
                write (6, '(A)') 'ERROR: init2_SAopt outcore and not give file'
                stop
            end if
            spin_f = 'SPINFILE'
        end if

        i0 = idx
        ! cycle for spin case
        do i = 1, n_spin
            if (spin_case(i) == 'AAAA') then
                call self%AAAA%init(iCore, poppa, popqa, popra, popsa, &
                                    ioppq, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(1, 1) = T
                if (self%incore) then
                    i0 = i0 + self%AAAA%len
                end if
            else if (spin_case(i) == 'AABB') then
                call self%AABB%init(iCore, poppa, popqa, poprb, popsb, &
                                    ioppq, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(1, 2) = T
                if (self%incore) then
                    i0 = i0 + self%AABB%len
                end if
            else if (spin_case(i) == 'AAAB') then
                call self%AAAB%init(iCore, poppa, popqa, popra, popsb, &
                                    ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(1, 3) = T
                if (self%incore) then
                    i0 = i0 + self%AAAB%len
                end if
            else if (spin_case(i) == 'AABA') then
                call self%AABA%init(iCore, poppa, popqa, poprb, popsa, &
                                    ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(1, 4) = T
                if (self%incore) then
                    i0 = i0 + self%AABA%len
                end if
            else if (spin_case(i) == 'BBAA') then
                call self%BBAA%init(iCore, poppb, popqb, popra, popsa, &
                                    ioppq, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(2, 1) = T
                if (self%incore) then
                    i0 = i0 + self%BBAA%len
                end if
            else if (spin_case(i) == 'BBBB') then
                call self%BBBB%init(iCore, poppb, popqb, poprb, popsb, &
                                    ioppq, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(2, 2) = T
                if (self%incore) then
                    i0 = i0 + self%BBBB%len
                end if
            else if (spin_case(i) == 'BBAB') then
                call self%BBAB%init(iCore, poppb, popqb, popra, popsb, &
                                    ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(2, 3) = T
                if (self%incore) then
                    i0 = i0 + self%BBAB%len
                end if
            else if (spin_case(i) == 'BBBA') then
                call self%BBBA%init(iCore, poppb, popqb, poprb, popsa, &
                                    ioppq, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(2, 4) = T
                if (self%incore) then
                    i0 = i0 + self%BBBA%len
                end if
            else if (spin_case(i) == 'ABAA') then
                call self%ABAA%init(iCore, poppa, popqb, popra, popsa, &
                                    1, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(3, 1) = T
                if (self%incore) then
                    i0 = i0 + self%ABAA%len
                end if
            else if (spin_case(i) == 'ABBB') then
                call self%ABBB%init(iCore, poppa, popqb, poprb, popsb, &
                                    1, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(3, 2) = T
                if (self%incore) then
                    i0 = i0 + self%ABBB%len
                end if
            else if (spin_case(i) == 'ABAB') then
                call self%ABAB%init(iCore, poppa, popqb, popra, popsb, &
                                    1, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(3, 3) = T
                if (self%incore) then
                    i0 = i0 + self%ABAB%len
                end if
            else if (spin_case(i) == 'ABBA') then
                call self%ABBA%init(iCore, poppa, popqb, poprb, popsa, &
                                    1, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(3, 4) = T
                if (self%incore) then
                    i0 = i0 + self%ABBA%len
                end if
            else if (spin_case(i) == 'BAAA') then
                call self%BAAA%init(iCore, poppb, popqa, popra, popsa, &
                                    1, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(4, 1) = T
                if (self%incore) then
                    i0 = i0 + self%BAAA%len
                end if
            else if (spin_case(i) == 'BABB') then
                call self%BABB%init(iCore, poppb, popqa, poprb, popsb, &
                                    1, ioprs, sym, self%incore, i0, spin_f(i))
                self%double_exist(4, 2) = T
                if (self%incore) then
                    i0 = i0 + self%BABB%len
                end if
            else if (spin_case(i) == 'BAAB') then
                call self%BAAB%init(iCore, poppb, popqa, popra, popsb, &
                                    1, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(4, 3) = T
                if (self%incore) then
                    i0 = i0 + self%BAAB%len
                end if
            else if (spin_case(i) == 'BABA') then
                call self%BABA%init(iCore, poppb, popqa, poprb, popsa, &
                                    1, 1, sym, self%incore, i0, spin_f(i))
                self%double_exist(4, 4) = T
                if (self%incore) then
                    i0 = i0 + self%BABA%len
                end if
            else
                write (6, '(2A)') 'ERROR: spin case in init2_SAopt error.', spin_case(i)
                stop
            end if
        end do
        self%size = i0 - idx

    end subroutine init2_SAopt

    ! 2027.7.3
    !> determine whether two population be the same
    function pop_eq(pop1, pop2) result(eq)
        implicit none

        integer*8, dimension(8), intent(in) :: pop1, pop2
        logical :: eq

        integer*8 :: i

        eq = .true.
        do i = 1, 8
            if (pop1(i) /= pop2(i)) then
                eq = .false.
                return
            end if
        end do

    end function pop_eq

    subroutine norm_SAopt(self, head, iCore, maxCore)
        implicit none

        class(SAopt), intent(in) :: self
        character(len=*), intent(in) :: head !> print as head//spin
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore

        integer*8 :: i, j

        if (self%single_exist(1)) then
            write (6, '(2A,F)') head, ' AA:', self%AA%norm(iCore, maxCore)
        end if
        if (self%single_exist(2)) then
            write (6, '(2A,F)') head, ' BB:', self%BB%norm(iCore, maxCore)
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                write (6, '(2A,F)') head, ' AB:', self%AB%norm(iCore, maxCore)
            end if
            if (self%single_exist(4)) then
                write (6, '(2A,F)') head, ' BA:', self%BA%norm(iCore, maxCore)
            end if
        end if
        if (self%double_exist(1, 1)) then
            write (6, '(2A,F)') head, ' AAAA:', self%AAAA%norm(iCore, maxCore)
        end if
        if (self%double_exist(1, 2)) then
            write (6, '(2A,F)') head, ' AABB:', self%AABB%norm(iCore, maxCore)
        end if
        if (self%SOC) then
            if (self%double_exist(1, 3)) then
                write (6, '(2A,F)') head, ' AAAB:', self%AAAB%norm(iCore, maxCore)
            end if
            if (self%double_exist(1, 4)) then
                write (6, '(2A,F)') head, ' AABA:', self%AABA%norm(iCore, maxCore)
            end if
        end if
        if (self%double_exist(2, 1)) then
            write (6, '(2A,F)') head, ' BBAA:', self%BBAA%norm(iCore, maxCore)
        end if
        if (self%double_exist(2, 2)) then
            write (6, '(2A,F)') head, ' BBBB:', self%BBBB%norm(iCore, maxCore)
        end if
        if (self%SOC) then
            if (self%double_exist(2, 3)) then
                write (6, '(2A,F)') head, ' BBAB:', self%BBAB%norm(iCore, maxCore)
            end if
            if (self%double_exist(2, 4)) then
                write (6, '(2A,F)') head, ' BBBA:', self%BBBA%norm(iCore, maxCore)
            end if
            if (self%double_exist(3, 1)) then
                write (6, '(2A,F)') head, ' ABAA:', self%ABAA%norm(iCore, maxCore)
            end if
            if (self%double_exist(3, 2)) then
                write (6, '(2A,F)') head, ' ABBB:', self%ABBB%norm(iCore, maxCore)
            end if
        end if
        if (self%double_exist(3, 3)) then
            write (6, '(2A,F)') head, ' ABAB:', self%ABAB%norm(iCore, maxCore)
        end if
        if (self%double_exist(3, 4)) then
            write (6, '(2A,F)') head, ' ABBA:', self%ABBA%norm(iCore, maxCore)
        end if
        if (self%SOC) then
            if (self%double_exist(4, 1)) then
                write (6, '(2A,F)') head, ' BAAA:', self%BAAA%norm(iCore, maxCore)
            end if
            if (self%double_exist(4, 2)) then
                write (6, '(2A,F)') head, ' BABB:', self%BABB%norm(iCore, maxCore)
            end if
        end if
        if (self%double_exist(4, 3)) then
            write (6, '(2A,F)') head, ' BAAB:', self%BAAB%norm(iCore, maxCore)
        end if
        if (self%double_exist(4, 4)) then
            write (6, '(2A,F)') head, ' BABA:', self%BABA%norm(iCore, maxCore)
        end if

    end subroutine norm_SAopt

    subroutine print_SAopt(self)
        implicit none

        class(SAopt), intent(in) :: self

        write (6, '(A)') 'irrep population:'
        write (6, '(A,8I5)') 'pa:', self%poppa
        write (6, '(A,8I5)') 'pb:', self%poppb
        write (6, '(A,8I5)') 'qa:', self%popqa
        write (6, '(A,8I5)') 'qb:', self%popqb
        write (6, '(A,8I5)') 'ra:', self%popra
        write (6, '(A,8I5)') 'rb:', self%poprb
        write (6, '(A,8I5)') 'sa:', self%popsa
        write (6, '(A,8I5)') 'sb:', self%popsb

        write (6, '(A,I2)') 'symmetry:', self%sym
        write (6, '(A,2I2)') 'ioppq/rs:', self%ioppq, self%ioprs

        write (6, '(A)') 'spin case:'
        if (self%single_exist(1)) then
            write (6, '(A)') 'AA'
        end if
        if (self%single_exist(2)) then
            write (6, '(A)') 'BB'
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                write (6, '(A)') 'AB'
            end if
            if (self%single_exist(4)) then
                write (6, '(A)') 'BA'
            end if
        end if
        if (self%double_exist(1, 1)) then
            write (6, '(A)') 'AAAA'
        end if
        if (self%double_exist(1, 2)) then
            write (6, '(A)') 'AABB'
        end if
        if (self%SOC) then
            if (self%double_exist(1, 3)) then
                write (6, '(A)') 'AAAB'
            end if
            if (self%double_exist(1, 4)) then
                write (6, '(A)') 'AABA'
            end if
        end if
        if (self%double_exist(2, 1)) then
            write (6, '(A)') 'BBAA'
        end if
        if (self%double_exist(2, 2)) then
            write (6, '(A)') 'BBBB'
        end if
        if (self%SOC) then
            if (self%double_exist(2, 3)) then
                write (6, '(A)') 'BBAB'
            end if
            if (self%double_exist(2, 4)) then
                write (6, '(A)') 'BBBA'
            end if
            if (self%double_exist(3, 1)) then
                write (6, '(A)') 'ABAA'
            end if
            if (self%double_exist(3, 2)) then
                write (6, '(A)') 'ABBB'
            end if
        end if
        if (self%double_exist(3, 3)) then
            write (6, '(A)') 'ABAB'
        end if
        if (self%double_exist(3, 4)) then
            write (6, '(A)') 'ABBA'
        end if
        if (self%SOC) then
            if (self%double_exist(4, 1)) then
                write (6, '(A)') 'BAAA'
            end if
            if (self%double_exist(4, 2)) then
                write (6, '(A)') 'BABB'
            end if
        end if
        if (self%double_exist(4, 3)) then
            write (6, '(A)') 'BAAB'
        end if
        if (self%double_exist(4, 4)) then
            write (6, '(A)') 'BABA'
        end if

        if (self%incore) then
            write (6, '(A,I15,F12.7,A)') 'size in memory:', self%size, &
                (self%size*8/1024.0d0**3), ' GiB'
        end if

    end subroutine print_SAopt

    function spin_exist_SAopt(self, spin_case) result(exist)
        implicit none

        class(SAopt), intent(in) :: self
        character(len=*), intent(in) :: spin_case
        logical :: exist

        exist = F
        select case (spin_case)
        case ('AA')
            if (self%single_exist(1)) then
                exist = T
            end if
        case ('BB')
            if (self%single_exist(2)) then
                exist = T
            end if
        case ('AAAA')
            if (self%double_exist(1, 1)) then
                exist = T
            end if
        case ('AABB')
            if (self%double_exist(1, 2)) then
                exist = T
            end if
        case ('BBAA')
            if (self%double_exist(2, 1)) then
                exist = T
            end if
        case ('BBBB')
            if (self%double_exist(2, 2)) then
                exist = T
            end if
        case ('ABAB')
            if (self%double_exist(3, 3)) then
                exist = T
            end if
        case ('ABBA')
            if (self%double_exist(3, 4)) then
                exist = T
            end if
        case ('BAAB')
            if (self%double_exist(4, 3)) then
                exist = T
            end if
        case ('BABA')
            if (self%double_exist(4, 4)) then
                exist = T
            end if
        end select
        if (self%SOC) then
            select case (spin_case)
            case ('AB')
                if (self%single_exist(3)) then
                    exist = T
                end if
            case ('BA')
                if (self%single_exist(4)) then
                    exist = T
                end if
            case ('AAAB')
                if (self%double_exist(1, 3)) then
                    exist = T
                end if
            case ('AABA')
                if (self%double_exist(1, 4)) then
                    exist = T
                end if
            case ('BBAB')
                if (self%double_exist(2, 3)) then
                    exist = T
                end if
            case ('BBBA')
                if (self%double_exist(2, 4)) then
                    exist = T
                end if
            case ('ABAA')
                if (self%double_exist(3, 1)) then
                    exist = T
                end if
            case ('ABBB')
                if (self%double_exist(3, 2)) then
                    exist = T
                end if
            case ('BAAA')
                if (self%double_exist(4, 1)) then
                    exist = T
                end if
            case ('BABB')
                if (self%double_exist(4, 2)) then
                    exist = T
                end if
            end select
        end if

        if (.not. exist) then
            write (6, '(2A)') 'WARN: require spin case not exist. ', spin_case
        end if

    end function spin_exist_SAopt

    subroutine arr_copy_SAopt(self, out)
        implicit none

        class(SAopt), intent(in) :: self
        type(SAopt), intent(inout) :: out

        if (self%single_exist(1) .and. out%single_exist(1)) then
            call out%AA%arr_copy(self%AA)
        end if
        if (self%single_exist(2) .and. out%single_exist(2)) then
            call out%BB%arr_copy(self%BB)
        end if
        if (self%double_exist(1, 1) .and. out%double_exist(1, 1)) then
            call out%AAAA%arr_copy(self%AAAA)
        end if
        if (self%double_exist(1, 2) .and. out%double_exist(1, 2)) then
            call out%AABB%arr_copy(self%AABB)
        end if
        if (self%double_exist(2, 1) .and. out%double_exist(2, 1)) then
            call out%BBAA%arr_copy(self%BBAA)
        end if
        if (self%double_exist(2, 2) .and. out%double_exist(2, 2)) then
            call out%BBBB%arr_copy(self%BBBB)
        end if
        if (self%double_exist(3, 3) .and. out%double_exist(3, 3)) then
            call out%ABAB%arr_copy(self%ABAB)
        end if
        if (self%double_exist(3, 4) .and. out%double_exist(3, 4)) then
            call out%ABBA%arr_copy(self%ABBA)
        end if
        if (self%double_exist(4, 3) .and. out%double_exist(4, 3)) then
            call out%BAAB%arr_copy(self%BAAB)
        end if
        if (self%double_exist(4, 4) .and. out%double_exist(4, 4)) then
            call out%BABA%arr_copy(self%BABA)
        end if
        if (self%SOC) then
            if (self%single_exist(3) .and. out%single_exist(3)) then
                call out%AB%arr_copy(self%AB)
            end if
            if (self%single_exist(4) .and. out%single_exist(4)) then
                call out%BA%arr_copy(self%BA)
            end if
            if (self%double_exist(1, 3) .and. out%double_exist(1, 3)) then
                call out%AAAB%arr_copy(self%AAAB)
            end if
            if (self%double_exist(1, 4) .and. out%double_exist(1, 4)) then
                call out%AABA%arr_copy(self%AABA)
            end if
            if (self%double_exist(2, 3) .and. out%double_exist(2, 3)) then
                call out%BBAB%arr_copy(self%BBAB)
            end if
            if (self%double_exist(2, 4) .and. out%double_exist(2, 4)) then
                call out%BBBA%arr_copy(self%BBBA)
            end if
            if (self%double_exist(3, 1) .and. out%double_exist(3, 1)) then
                call out%ABAA%arr_copy(self%ABAA)
            end if
            if (self%double_exist(3, 2) .and. out%double_exist(3, 2)) then
                call out%ABBB%arr_copy(self%ABBB)
            end if
            if (self%double_exist(4, 1) .and. out%double_exist(4, 1)) then
                call out%BAAA%arr_copy(self%BAAA)
            end if
            if (self%double_exist(4, 2) .and. out%double_exist(4, 2)) then
                call out%BABB%arr_copy(self%BABB)
            end if
        end if

    end subroutine arr_copy_SAopt

    subroutine arr_clear_SAopt(self)
        implicit none

        class(SAopt), intent(inout) :: self

        if (self%single_exist(1)) then
            call dscal(self%AA%len, 0.0d0, self%AA%arr, 1)
        end if
        if (self%single_exist(2)) then
            call dscal(self%BB%len, 0.0d0, self%BB%arr, 1)
        end if
        if (self%double_exist(1, 1)) then
            call dscal(self%AAAA%len, 0.0d0, self%AAAA%arr, 1)
        end if
        if (self%double_exist(1, 2)) then
            call dscal(self%AABB%len, 0.0d0, self%AABB%arr, 1)
        end if
        if (self%double_exist(2, 1)) then
            call dscal(self%BBAA%len, 0.0d0, self%BBAA%arr, 1)
        end if
        if (self%double_exist(2, 2)) then
            call dscal(self%BBBB%len, 0.0d0, self%BBBB%arr, 1)
        end if
        if (self%double_exist(3, 3)) then
            call dscal(self%ABAB%len, 0.0d0, self%ABAB%arr, 1)
        end if
        if (self%double_exist(3, 4)) then
            call dscal(self%ABBA%len, 0.0d0, self%ABBA%arr, 1)
        end if
        if (self%double_exist(4, 3)) then
            call dscal(self%BAAB%len, 0.0d0, self%BAAB%arr, 1)
        end if
        if (self%double_exist(4, 4)) then
            call dscal(self%BABA%len, 0.0d0, self%BABA%arr, 1)
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                call dscal(self%AB%len, 0.0d0, self%AB%arr, 1)
            end if
            if (self%single_exist(4)) then
                call dscal(self%BA%len, 0.0d0, self%BA%arr, 1)
            end if
            if (self%double_exist(1, 3)) then
                call dscal(self%AAAB%len, 0.0d0, self%AAAB%arr, 1)
            end if
            if (self%double_exist(1, 4)) then
                call dscal(self%AABA%len, 0.0d0, self%AABA%arr, 1)
            end if
            if (self%double_exist(2, 3)) then
                call dscal(self%BBAB%len, 0.0d0, self%BBAB%arr, 1)
            end if
            if (self%double_exist(2, 4)) then
                call dscal(self%BBBA%len, 0.0d0, self%BBBA%arr, 1)
            end if
            if (self%double_exist(3, 1)) then
                call dscal(self%ABAA%len, 0.0d0, self%ABAA%arr, 1)
            end if
            if (self%double_exist(3, 2)) then
                call dscal(self%ABBB%len, 0.0d0, self%ABBB%arr, 1)
            end if
            if (self%double_exist(4, 1)) then
                call dscal(self%BAAA%len, 0.0d0, self%BAAA%arr, 1)
            end if
            if (self%double_exist(4, 2)) then
                call dscal(self%BABB%len, 0.0d0, self%BABB%arr, 1)
            end if
        end if

    end subroutine arr_clear_SAopt

    subroutine copy_SAopt(self, out, iCore, idx)
        implicit none

        class(SAopt), intent(in) :: self
        type(SAopt), intent(out) :: out
        real*8, dimension(*), intent(inout) :: iCore !> work space
        integer*8, intent(in) :: idx !> poistion of out

        integer*8 :: i0

        i0 = idx

        out%poppa = self%poppa
        out%poppb = self%poppb
        out%popqa = self%popqa
        out%popqb = self%popqb
        out%ioppq = self%ioppq
        out%popra = self%popra
        out%poprb = self%poprb
        out%popsa = self%popsa
        out%popsb = self%popsb
        out%ioprs = self%ioprs
        out%sym = self%sym
        out%UHF = self%UHF
        out%single_exist = self%single_exist
        out%double_exist = self%double_exist
        out%particle = self%particle
        out%incore = self%incore

        if (self%single_exist(1)) then
            call copy_opt(self%AA, out%AA, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%AA%len
            end if
        end if
        if (self%single_exist(2)) then
            call copy_opt(self%BB, out%BB, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%BB%len
            end if
        end if
        if (self%double_exist(1, 1)) then
            call copy_opt(self%AAAA, out%AAAA, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%AAAA%len
            end if
        end if
        if (self%double_exist(1, 2)) then
            call copy_opt(self%AABB, out%AABB, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%AABB%len
            end if
        end if
        if (self%double_exist(2, 1)) then
            call copy_opt(self%BBAA, out%BBAA, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%BBAA%len
            end if
        end if
        if (self%double_exist(2, 2)) then
            call copy_opt(self%BBBB, out%BBBB, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%BBBB%len
            end if
        end if
        if (self%double_exist(3, 3)) then
            call copy_opt(self%ABAB, out%ABAB, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%ABAB%len
            end if
        end if
        if (self%double_exist(3, 4)) then
            call copy_opt(self%ABBA, out%ABBA, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%ABBA%len
            end if
        end if
        if (self%double_exist(4, 3)) then
            call copy_opt(self%BAAB, out%BAAB, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%BAAB%len
            end if
        end if
        if (self%double_exist(4, 4)) then
            call copy_opt(self%BABA, out%BABA, out%incore, iCore, i0)
            if (out%incore) then
                i0 = i0 + out%BABA%len
            end if
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                call copy_opt(self%AB, out%AB, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%AB%len
                end if
            end if
            if (self%single_exist(4)) then
                call copy_opt(self%BA, out%BA, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%BA%len
                end if
            end if
            if (self%double_exist(1, 3)) then
                call copy_opt(self%AAAB, out%AAAB, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%AAAB%len
                end if
            end if
            if (self%double_exist(1, 4)) then
                call copy_opt(self%AABA, out%AABA, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%AABA%len
                end if
            end if
            if (self%double_exist(2, 3)) then
                call copy_opt(self%BBAB, out%BBAB, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%BBAB%len
                end if
            end if
            if (self%double_exist(2, 4)) then
                call copy_opt(self%BBBA, out%BBBA, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%BBBA%len
                end if
            end if
            if (self%double_exist(3, 1)) then
                call copy_opt(self%ABAA, out%ABAA, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%ABAA%len
                end if
            end if
            if (self%double_exist(3, 2)) then
                call copy_opt(self%ABBB, out%ABBB, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%ABBB%len
                end if
            end if
            if (self%double_exist(4, 1)) then
                call copy_opt(self%BAAA, out%BAAA, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%BAAA%len
                end if
            end if
            if (self%double_exist(4, 2)) then
                call copy_opt(self%BABB, out%BABB, out%incore, iCore, i0)
                if (out%incore) then
                    i0 = i0 + out%BABB%len
                end if
            end if
        end if

        out%size = i0 - idx
        if (self%incore) then
            if (out%size /= self%size) then
                write (6, '(A)') 'ERROR: size not correct in copy_SAopt.'
                stop
            end if
        end if

    end subroutine copy_SAopt

    subroutine read_SAopt(self)
        implicit none

        class(SAopt), intent(in) :: self

        if (self%single_exist(1)) then
            call getallso(self%AA%arr, 1, self%AA%f)
        end if
        if (self%single_exist(2)) then
            call getallso(self%BB%arr, 1, self%BB%f)
        end if
        if (self%double_exist(1, 1)) then
            call getallso(self%AAAA%arr, 1, self%AAAA%f)
        end if
        if (self%double_exist(1, 2)) then
            call getallso(self%AABB%arr, 1, self%AABB%f)
        end if
        if (self%double_exist(2, 1)) then
            call getallso(self%BBAA%arr, 1, self%BBAA%f)
        end if
        if (self%double_exist(2, 2)) then
            call getallso(self%BBBB%arr, 1, self%BBBB%f)
        end if
        if (self%double_exist(3, 3)) then
            call getallso(self%ABAB%arr, 1, self%ABAB%f)
        end if
        if (self%double_exist(3, 4)) then
            call getallso(self%ABBA%arr, 1, self%ABBA%f)
        end if
        if (self%double_exist(4, 3)) then
            call getallso(self%BAAB%arr, 1, self%BAAB%f)
        end if
        if (self%double_exist(4, 4)) then
            call getallso(self%BABA%arr, 1, self%BABA%f)
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                call getallso(self%AB%arr, 1, self%AB%f)
            end if
            if (self%single_exist(4)) then
                call getallso(self%BA%arr, 1, self%BA%f)
            end if
            if (self%double_exist(1, 3)) then
                call getallso(self%AAAB%arr, 1, self%AAAB%f)
            end if
            if (self%double_exist(1, 4)) then
                call getallso(self%AABA%arr, 1, self%AABA%f)
            end if
            if (self%double_exist(2, 3)) then
                call getallso(self%BBAB%arr, 1, self%BBAB%f)
            end if
            if (self%double_exist(2, 4)) then
                call getallso(self%BBBA%arr, 1, self%BBBA%f)
            end if
            if (self%double_exist(3, 1)) then
                call getallso(self%ABAA%arr, 1, self%ABAA%f)
            end if
            if (self%double_exist(3, 2)) then
                call getallso(self%ABBB%arr, 1, self%ABBB%f)
            end if
            if (self%double_exist(4, 1)) then
                call getallso(self%BAAA%arr, 1, self%BAAA%f)
            end if
            if (self%double_exist(4, 2)) then
                call getallso(self%BABB%arr, 1, self%BABB%f)
            end if
        end if

    end subroutine read_SAopt

    subroutine write_SAopt(self)
        implicit none

        class(SAopt), intent(in) :: self

        if (self%single_exist(1)) then
            call putallso(self%AA%arr, 1, self%AA%f)
        end if
        if (self%single_exist(2)) then
            call putallso(self%BB%arr, 1, self%BB%f)
        end if
        if (self%double_exist(1, 1)) then
            call putallso(self%AAAA%arr, 1, self%AAAA%f)
        end if
        if (self%double_exist(1, 2)) then
            call putallso(self%AABB%arr, 1, self%AABB%f)
        end if
        if (self%double_exist(2, 1)) then
            call putallso(self%BBAA%arr, 1, self%BBAA%f)
        end if
        if (self%double_exist(2, 2)) then
            call putallso(self%BBBB%arr, 1, self%BBBB%f)
        end if
        if (self%double_exist(3, 3)) then
            call putallso(self%ABAB%arr, 1, self%ABAB%f)
        end if
        if (self%double_exist(3, 4)) then
            call putallso(self%ABBA%arr, 1, self%ABBA%f)
        end if
        if (self%double_exist(4, 3)) then
            call putallso(self%BAAB%arr, 1, self%BAAB%f)
        end if
        if (self%double_exist(4, 4)) then
            call putallso(self%BABA%arr, 1, self%BABA%f)
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                call putallso(self%AB%arr, 1, self%AB%f)
            end if
            if (self%single_exist(4)) then
                call putallso(self%BA%arr, 1, self%BA%f)
            end if
            if (self%double_exist(1, 3)) then
                call putallso(self%AAAB%arr, 1, self%AAAB%f)
            end if
            if (self%double_exist(1, 4)) then
                call putallso(self%AABA%arr, 1, self%AABA%f)
            end if
            if (self%double_exist(2, 3)) then
                call putallso(self%BBAB%arr, 1, self%BBAB%f)
            end if
            if (self%double_exist(2, 4)) then
                call putallso(self%BBBA%arr, 1, self%BBBA%f)
            end if
            if (self%double_exist(3, 1)) then
                call putallso(self%ABAA%arr, 1, self%ABAA%f)
            end if
            if (self%double_exist(3, 2)) then
                call putallso(self%ABBB%arr, 1, self%ABBB%f)
            end if
            if (self%double_exist(4, 1)) then
                call putallso(self%BAAA%arr, 1, self%BAAA%f)
            end if
            if (self%double_exist(4, 2)) then
                call putallso(self%BABB%arr, 1, self%BABB%f)
            end if
        end if

    end subroutine write_SAopt

    subroutine create_f_SAopt(self, iCore)
        implicit none

        class(SAopt), intent(in) :: self
        real*8, dimension(*), intent(inout) :: iCore !> free space

        if (self%single_exist(1)) then
            call clmoioso0(self%AA%sym, iCore, self%AA%popp, self%AA%popq, self%AA%ioppq, &
                           self%AA%popr, self%AA%pops, self%AA%ioprs, self%AA%f)
        end if
        if (self%single_exist(2)) then
            call clmoioso0(self%BB%sym, iCore, self%BB%popp, self%BB%popq, self%BB%ioppq, &
                           self%BB%popr, self%BB%pops, self%BB%ioprs, self%BB%f)
        end if
        if (self%double_exist(1, 1)) then
            call clmoioso0(self%AAAA%sym, iCore, self%AAAA%popp, self%AAAA%popq, self%AAAA%ioppq, &
                           self%AAAA%popr, self%AAAA%pops, self%AAAA%ioprs, self%AAAA%f)
        end if
        if (self%double_exist(1, 2)) then
            call clmoioso0(self%AABB%sym, iCore, self%AABB%popp, self%AABB%popq, self%AABB%ioppq, &
                           self%AABB%popr, self%AABB%pops, self%AABB%ioprs, self%AABB%f)
        end if
        if (self%double_exist(2, 1)) then
            call clmoioso0(self%BBAA%sym, iCore, self%BBAA%popp, self%BBAA%popq, self%BBAA%ioppq, &
                           self%BBAA%popr, self%BBAA%pops, self%BBAA%ioprs, self%BBAA%f)
        end if
        if (self%double_exist(2, 2)) then
            call clmoioso0(self%BBBB%sym, iCore, self%BBBB%popp, self%BBBB%popq, self%BBBB%ioppq, &
                           self%BBBB%popr, self%BBBB%pops, self%BBBB%ioprs, self%BBBB%f)
        end if
        if (self%double_exist(3, 3)) then
            call clmoioso0(self%ABAB%sym, iCore, self%ABAB%popp, self%ABAB%popq, self%ABAB%ioppq, &
                           self%ABAB%popr, self%ABAB%pops, self%ABAB%ioprs, self%ABAB%f)
        end if
        if (self%double_exist(3, 4)) then
            call clmoioso0(self%ABBA%sym, iCore, self%ABBA%popp, self%ABBA%popq, self%ABBA%ioppq, &
                           self%ABBA%popr, self%ABBA%pops, self%ABBA%ioprs, self%ABBA%f)
        end if
        if (self%double_exist(4, 3)) then
            call clmoioso0(self%BAAB%sym, iCore, self%BAAB%popp, self%BAAB%popq, self%BAAB%ioppq, &
                           self%BAAB%popr, self%BAAB%pops, self%BAAB%ioprs, self%BAAB%f)
        end if
        if (self%double_exist(4, 4)) then
            call clmoioso0(self%BABA%sym, iCore, self%BABA%popp, self%BABA%popq, self%BABA%ioppq, &
                           self%BABA%popr, self%BABA%pops, self%BABA%ioprs, self%BABA%f)
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                call clmoioso0(self%AB%sym, iCore, self%AB%popp, self%AB%popq, self%AB%ioppq, &
                               self%AB%popr, self%AB%pops, self%AB%ioprs, self%AB%f)
            end if
            if (self%single_exist(4)) then
                call clmoioso0(self%BA%sym, iCore, self%BA%popp, self%BA%popq, self%BA%ioppq, &
                               self%BA%popr, self%BA%pops, self%BA%ioprs, self%BA%f)
            end if
            if (self%double_exist(1, 3)) then
                call clmoioso0(self%AAAB%sym, iCore, self%AAAB%popp, self%AAAB%popq, self%AAAB%ioppq, &
                               self%AAAB%popr, self%AAAB%pops, self%AAAB%ioprs, self%AAAB%f)
            end if
            if (self%double_exist(1, 4)) then
                call clmoioso0(self%AABA%sym, iCore, self%AABA%popp, self%AABA%popq, self%AABA%ioppq, &
                               self%AABA%popr, self%AABA%pops, self%AABA%ioprs, self%AABA%f)
            end if
            if (self%double_exist(2, 3)) then
                call clmoioso0(self%BBAB%sym, iCore, self%BBAB%popp, self%BBAB%popq, self%BBAB%ioppq, &
                               self%BBAB%popr, self%BBAB%pops, self%BBAB%ioprs, self%BBAB%f)
            end if
            if (self%double_exist(2, 4)) then
                call clmoioso0(self%BBBA%sym, iCore, self%BBBA%popp, self%BBBA%popq, self%BBBA%ioppq, &
                               self%BBBA%popr, self%BBBA%pops, self%BBBA%ioprs, self%BBBA%f)
            end if
            if (self%double_exist(3, 1)) then
                call clmoioso0(self%ABAA%sym, iCore, self%ABAA%popp, self%ABAA%popq, self%ABAA%ioppq, &
                               self%ABAA%popr, self%ABAA%pops, self%ABAA%ioprs, self%ABAA%f)
            end if
            if (self%double_exist(3, 2)) then
                call clmoioso0(self%ABBB%sym, iCore, self%ABBB%popp, self%ABBB%popq, self%ABBB%ioppq, &
                               self%ABBB%popr, self%ABBB%pops, self%ABBB%ioprs, self%ABBB%f)
            end if
            if (self%double_exist(4, 1)) then
                call clmoioso0(self%BAAA%sym, iCore, self%BAAA%popp, self%BAAA%popq, self%BAAA%ioppq, &
                               self%BAAA%popr, self%BAAA%pops, self%BAAA%ioprs, self%BAAA%f)
            end if
            if (self%double_exist(4, 2)) then
                call clmoioso0(self%BABB%sym, iCore, self%BABB%popp, self%BABB%popq, self%BABB%ioppq, &
                               self%BABB%popr, self%BABB%pops, self%BABB%ioprs, self%BABB%f)
            end if
        end if

    end subroutine create_f_SAopt

    subroutine add_SAopt(self, in, fact)
        implicit none

        class(SAopt), intent(inout) :: self
        class(SAopt), intent(in) :: in
        real*8, intent(in) :: fact

        if (self%single_exist(1)) then
            call daxpy(self%AA%len, fact, in%AA%arr, 1, self%AA%arr, 1)
        end if
        if (self%single_exist(2)) then
            call daxpy(self%BB%len, fact, in%BB%arr, 1, self%BB%arr, 1)
        end if
        if (self%double_exist(1, 1)) then
            call daxpy(self%AAAA%len, fact, in%AAAA%arr, 1, self%AAAA%arr, 1)
        end if
        if (self%double_exist(1, 2)) then
            call daxpy(self%AABB%len, fact, in%AABB%arr, 1, self%AABB%arr, 1)
        end if
        if (self%double_exist(2, 1)) then
            call daxpy(self%BBAA%len, fact, in%BBAA%arr, 1, self%BBAA%arr, 1)
        end if
        if (self%double_exist(2, 2)) then
            call daxpy(self%BBBB%len, fact, in%BBBB%arr, 1, self%BBBB%arr, 1)
        end if
        if (self%double_exist(3, 3)) then
            call daxpy(self%ABAB%len, fact, in%ABAB%arr, 1, self%ABAB%arr, 1)
        end if
        if (self%double_exist(3, 4)) then
            call daxpy(self%ABBA%len, fact, in%ABBA%arr, 1, self%ABBA%arr, 1)
        end if
        if (self%double_exist(4, 3)) then
            call daxpy(self%BAAB%len, fact, in%BAAB%arr, 1, self%BAAB%arr, 1)
        end if
        if (self%double_exist(4, 4)) then
            call daxpy(self%BABA%len, fact, in%BABA%arr, 1, self%BABA%arr, 1)
        end if
        if (self%SOC) then
            if (self%single_exist(3)) then
                call daxpy(self%AB%len, fact, in%AB%arr, 1, self%AB%arr, 1)
            end if
            if (self%single_exist(4)) then
                call daxpy(self%BA%len, fact, in%BA%arr, 1, self%BA%arr, 1)
            end if
            if (self%double_exist(1, 3)) then
                call daxpy(self%AAAB%len, fact, in%AAAB%arr, 1, self%AAAB%arr, 1)
            end if
            if (self%double_exist(1, 4)) then
                call daxpy(self%AABA%len, fact, in%AABA%arr, 1, self%AABA%arr, 1)
            end if
            if (self%double_exist(2, 3)) then
                call daxpy(self%BBAB%len, fact, in%BBAB%arr, 1, self%BBAB%arr, 1)
            end if
            if (self%double_exist(2, 4)) then
                call daxpy(self%BBBA%len, fact, in%BBBA%arr, 1, self%BBBA%arr, 1)
            end if
            if (self%double_exist(3, 1)) then
                call daxpy(self%ABAA%len, fact, in%ABAA%arr, 1, self%ABAA%arr, 1)
            end if
            if (self%double_exist(3, 2)) then
                call daxpy(self%ABBB%len, fact, in%ABBB%arr, 1, self%ABBB%arr, 1)
            end if
            if (self%double_exist(4, 1)) then
                call daxpy(self%BAAA%len, fact, in%BAAA%arr, 1, self%BAAA%arr, 1)
            end if
            if (self%double_exist(4, 2)) then
                call daxpy(self%BABB%len, fact, in%BABB%arr, 1, self%BABB%arr, 1)
            end if
        end if

    end subroutine add_SAopt

    subroutine transpose_SAopt(self, out)
        implicit none

        class(SAopt), intent(in) :: self
        class(SAopt), intent(inout) :: out

        if (self%double_exist(1, 1)) then
            call itranspso(self%AAAA%arr, out%AAAA%arr, self%AAAA%popp, self%AAAA%popq, self%AAAA%ioppq, &
                           self%AAAA%popr, self%AAAA%pops, self%AAAA%ioprs, self%AAAA%sym)
        end if
        if (self%double_exist(1, 2)) then
            call itranspso(self%AABB%arr, out%BBAA%arr, self%AABB%popp, self%AABB%popq, self%AABB%ioppq, &
                           self%AABB%popr, self%AABB%pops, self%AABB%ioprs, self%AABB%sym)
        end if
        if (self%double_exist(2, 1)) then
            call itranspso(self%BBAA%arr, out%AABB%arr, self%BBAA%popp, self%BBAA%popq, self%BBAA%ioppq, &
                           self%BBAA%popr, self%BBAA%pops, self%BBAA%ioprs, self%BBAA%sym)
        end if
        if (self%double_exist(2, 2)) then
            call itranspso(self%BBBB%arr, out%BBBB%arr, self%BBBB%popp, self%BBBB%popq, self%BBBB%ioppq, &
                           self%BBBB%popr, self%BBBB%pops, self%BBBB%ioprs, self%BBBB%sym)
        end if
        if (self%double_exist(3, 3)) then
            call itranspso(self%ABAB%arr, out%ABAB%arr, self%ABAB%popp, self%ABAB%popq, self%ABAB%ioppq, &
                           self%ABAB%popr, self%ABAB%pops, self%ABAB%ioprs, self%ABAB%sym)
        end if
        if (self%double_exist(3, 4)) then
            call itranspso(self%ABBA%arr, out%BAAB%arr, self%ABBA%popp, self%ABBA%popq, self%ABBA%ioppq, &
                           self%ABBA%popr, self%ABBA%pops, self%ABBA%ioprs, self%ABBA%sym)
        end if
        if (self%double_exist(4, 3)) then
            call itranspso(self%BAAB%arr, out%ABBA%arr, self%BAAB%popp, self%BAAB%popq, self%BAAB%ioppq, &
                           self%BAAB%popr, self%BAAB%pops, self%BAAB%ioprs, self%BAAB%sym)
        end if
        if (self%double_exist(4, 4)) then
            call itranspso(self%BABA%arr, out%BABA%arr, self%BABA%popp, self%BABA%popq, self%BABA%ioppq, &
                           self%BABA%popr, self%BABA%pops, self%BABA%ioprs, self%BABA%sym)
        end if
        if (self%SOC) then
            if (self%double_exist(1, 3)) then
                call itranspso(self%AAAB%arr, out%ABAA%arr, self%AAAB%popp, self%AAAB%popq, self%AAAB%ioppq, &
                               self%AAAB%popr, self%AAAB%pops, self%AAAB%ioprs, self%AAAB%sym)
            end if
            if (self%double_exist(1, 4)) then
                call itranspso(self%AABA%arr, out%BAAA%arr, self%AABA%popp, self%AABA%popq, self%AABA%ioppq, &
                               self%AABA%popr, self%AABA%pops, self%AABA%ioprs, self%AABA%sym)
            end if
            if (self%double_exist(2, 3)) then
                call itranspso(self%BBAB%arr, out%ABBB%arr, self%BBAB%popp, self%BBAB%popq, self%BBAB%ioppq, &
                               self%BBAB%popr, self%BBAB%pops, self%BBAB%ioprs, self%BBAB%sym)
            end if
            if (self%double_exist(2, 4)) then
                call itranspso(self%BBBA%arr, out%BABB%arr, self%BBBA%popp, self%BBBA%popq, self%BBBA%ioppq, &
                               self%BBBA%popr, self%BBBA%pops, self%BBBA%ioprs, self%BBBA%sym)
            end if
            if (self%double_exist(3, 1)) then
                call itranspso(self%ABAA%arr, out%AAAB%arr, self%ABAA%popp, self%ABAA%popq, self%ABAA%ioppq, &
                               self%ABAA%popr, self%ABAA%pops, self%ABAA%ioprs, self%ABAA%sym)
            end if
            if (self%double_exist(3, 2)) then
                call itranspso(self%ABBB%arr, out%BBAB%arr, self%ABBB%popp, self%ABBB%popq, self%ABBB%ioppq, &
                               self%ABBB%popr, self%ABBB%pops, self%ABBB%ioprs, self%ABBB%sym)
            end if
            if (self%double_exist(4, 1)) then
                call itranspso(self%BAAA%arr, out%AABA%arr, self%BAAA%popp, self%BAAA%popq, self%BAAA%ioppq, &
                               self%BAAA%popr, self%BAAA%pops, self%BAAA%ioprs, self%BAAA%sym)
            end if
            if (self%double_exist(4, 2)) then
                call itranspso(self%BABB%arr, out%BBBA%arr, self%BABB%popp, self%BABB%popq, self%BABB%ioppq, &
                               self%BABB%popr, self%BABB%pops, self%BABB%ioprs, self%BABB%sym)
            end if
        end if

    end subroutine transpose_SAopt

end module mod_operator
