!> Such module is an interface to run spin adapted contracton
module mod_mult
    use mod_operator
    implicit none
    public

contains
    ! type=1: t(ab,ij) += fact*P(ab){t(a,i)*t(b,j)}
    ! type=2: t(ab,ij) += fact*P(ab|ij){t(a,i)*f(b,j)}
    ! type=3: t(ab,ij) += fact*P(ab){t(a,i)*f(b,j)}
    ! type=4: t(ab,ij) += fact*P(ij){t(a,i)*f(b,j)}
    subroutine effectT2_SA(T2, T1, fact, iCore, maxCore, type, R1)
        implicit none

        type(SAopt), intent(inout) :: T2
        type(SAopt), intent(in) :: T1
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        integer*8, intent(in), optional :: type !> P(ab)/P(ab|ij)/P(ab)/P(ij)
        type(SAopt), intent(in), optional :: R1

        logical :: spin_exist = T

        spin_exist = spin_exist .and. T2%spin_exist('AAAA')
        spin_exist = spin_exist .and. T2%spin_exist('ABAB')
        spin_exist = spin_exist .and. T1%spin_exist('AA')

        if ((.not. present(type)) .or. (type == 1)) then
            if (T2%UHF) then
                spin_exist = spin_exist .and. T2%spin_exist('BBBB')
                spin_exist = spin_exist .and. T1%spin_exist('BB')

                call effectT2_U(T2%AAAA, T2%ABAB, T2%BBBB, T2%ABAB, &
                                T1%AA, T1%BB, T1%AA, T1%BB, &
                                fact, iCore, maxCore, 1)
            else
                call effectT2_R(T2%AAAA, T2%ABAB, T1%AA, T1%AA, fact, iCore, maxCore, 1)
            end if
        else if (present(type) .and. (type /= 1)) then
            if (.not. present(R1)) then
                write (6, '(A)') 'ERROR: effectT2_SA with type/=1 must give a R1.'
                stop
            end if
            spin_exist = spin_exist .and. R1%spin_exist('AA')

            if (T2%UHF) then
                spin_exist = spin_exist .and. T2%spin_exist('BBBB')
                spin_exist = spin_exist .and. T1%spin_exist('BB')
                spin_exist = spin_exist .and. R1%spin_exist('BB')

                call effectT2_U(T2%AAAA, T2%ABAB, T2%BBBB, T2%BABA, &
                                T1%AA, T1%BB, R1%AA, R1%BB, &
                                fact, iCore, maxCore, type)
            else
                call effectT2_R(T2%AAAA, T2%ABAB, T1%AA, R1%AA, fact, iCore, maxCore, type)
            end if
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in effectT2_SA error.'
            stop
        end if

    end subroutine effectT2_SA

    !> c(a,i) += fact*r(e,i)*F(a,e)
    subroutine c1_r1Fae_SA(c1, r1, Fae, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1
        type(SAopt), intent(in) :: r1, Fae
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. Fae%spin_exist('AA')

        if (c1%UHF) then
            spin_exist = spin_exist .and. c1%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. Fae%spin_exist('BB')

            call c1_r1Fae_U(c1%AA, c1%BB, r1%AA, r1%BB, Fae%AA, Fae%BB, &
                            fact, iCore, maxCore)
        else
            call c1_r1Fae_R(c1%AA, r1%AA, Fae%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r1Fae_SA error.'
            stop
        end if

    end subroutine c1_r1Fae_SA

    !> c(a,i) -= fact*r(a,m)*F(m,i)
    subroutine c1_r1Fmi_SA(c1, r1, Fmi, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1
        type(SAopt), intent(in) :: r1, Fmi
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. Fmi%spin_exist('AA')

        if (c1%UHF) then
            spin_exist = spin_exist .and. c1%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. Fmi%spin_exist('BB')

            call c1_r1Fmi_U(c1%AA, c1%BB, r1%AA, r1%BB, Fmi%AA, Fmi%BB, &
                            fact, iCore, maxCore)
        else
            call c1_r1Fmi_R(c1%AA, r1%AA, Fmi%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r1Fmi_SA error.'
            stop
        end if

    end subroutine c1_r1Fmi_SA

    !> c(a,i) += fact*r(ae,im)*F(e,m)
    subroutine c1_r2Fme_SA(c1, r2, Fme, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1
        type(SAopt), intent(in) :: r2, Fme
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1%spin_exist('AA')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fme%spin_exist('AA')

        if (c1%UHF) then
            spin_exist = spin_exist .and. c1%spin_exist('BB')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fme%spin_exist('BB')

            call c1_r2Fme_U(c1%AA, c1%BB, r2%AAAA, r2%ABAB, r2%BBBB, Fme%AA, Fme%BB, &
                            fact, iCore, maxCore)
        else
            call c1_r2Fme_R(c1%AA, r2%AAAA, r2%ABAB, Fme%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r2Fme_SA error.'
            stop
        end if

    end subroutine c1_r2Fme_SA

    !> c(a,i) += fact*r(f,n)*Vnafi(fn,ai), need V(ia,bj)=V(jb,ai)
    subroutine c1_r1iabj_SA(c1, r1, iabj, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1
        type(SAopt), intent(in) :: r1, iabj
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. iabj%spin_exist('AAAA')
        spin_exist = spin_exist .and. iabj%spin_exist('AABB')

        if (c1%UHF) then
            spin_exist = spin_exist .and. c1%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. iabj%spin_exist('BBBB')

            call c1_r1iabj_U(c1%AA, c1%BB, r1%AA, r1%BB, iabj%AAAA, iabj%AABB, iabj%BBBB, &
                             fact, iCore, maxCore)
        else
            call c1_r1iabj_R(c1%AA, r1%AA, iabj%AAAA, iabj%AABB, &
                             fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r1iabj_SA error.'
            stop
        end if

    end subroutine c1_r1iabj_SA

    !> c(a,i) += fact*0.5*r(fe,im)*V(fe,am)
    subroutine c1_r2abci_SA(c1, r2, abci, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1
        type(SAopt), intent(in) :: r2, abci
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1%spin_exist('AA')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abci%spin_exist('AAAA')
        spin_exist = spin_exist .and. abci%spin_exist('ABAB')

        if (c1%UHF) then
            spin_exist = spin_exist .and. c1%spin_exist('BB')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BABA')

            call c1_r2abci_U(c1%AA, c1%BB, r2%AAAA, r2%ABAB, r2%BBBB, &
                             abci%AAAA, abci%ABAB, abci%BBBB, abci%BABA, &
                             fact, iCore, maxCore)
        else
            call c1_r2abci_R(c1%AA, r2%AAAA, r2%ABAB, abci%AAAA, abci%ABAB, &
                             fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r2abci_SA error.'
            stop
        end if

    end subroutine c1_r2abci_SA

    !> c(a,i) -= fact*0.5*r(ae,mn)*V(mn,ie)
    subroutine c1_r2ijka_SA(c1, r2, ijka, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1
        type(SAopt), intent(in) :: r2, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1%spin_exist('AA')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (c1%UHF) then
            spin_exist = spin_exist .and. c1%spin_exist('BB')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call c1_r2ijka_U(c1%AA, c1%BB, r2%AAAA, r2%ABAB, r2%BBBB, &
                             ijka%AAAA, ijka%ABAB, ijka%BBBB, ijka%BABA, &
                             fact, iCore, maxCore)
        else
            call c1_r2ijka_R(c1%AA, r2%AAAA, r2%ABAB, ijka%AAAA, ijka%ABAB, &
                             fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r2ijka error.'
            stop
        end if

    end subroutine c1_r2ijka_SA

    !> c(ab,ij) += fact*P(ab){r(ae,ij)*F(b,e)}
    subroutine c2_r2Fae_SA(c2, r2, Fae, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r2, Fae
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fae%spin_exist('AA')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fae%spin_exist('BB')

            call c2_r2Fae_U(c2%AAAA, c2%ABAB, c2%BBBB, r2%AAAA, r2%ABAB, r2%BBBB, &
                            Fae%AA, Fae%BB, fact, iCore, maxCore)
        else
            call c2_r2Fae_R(c2%ABAB, r2%ABAB, Fae%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1_r2Fae error.'
            stop
        end if

    end subroutine c2_r2Fae_SA

    !> c(ab,ij) -= fact*P(ij){r(ab,im)*F(m,j)}
    subroutine c2_r2Fmi_SA(c2, r2, Fmi, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r2, Fmi
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fmi%spin_exist('AA')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fmi%spin_exist('BB')

            call c2_r2Fmi_U(c2%AAAA, c2%ABAB, c2%BBBB, r2%AAAA, r2%ABAB, r2%BBBB, &
                            Fmi%AA, Fmi%BB, fact, iCore, maxCore)
        else
            call c2_r2Fmi_R(c2%ABAB, r2%ABAB, Fmi%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2_r2Fmi_SA error.'
            stop
        end if

    end subroutine c2_r2Fmi_SA

    !> c(ab,ij) += 0.5*fact*r(ab,mn)*W(mn,ij)
    subroutine c2_r2Wmnij_SA(c2, r2, Wmnij, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r2, Wmnij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Wmnij%spin_exist('ABAB')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wmnij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wmnij%spin_exist('BBBB')

            call c2_r2Wmnij_U(c2%AAAA, c2%ABAB, c2%BBBB, r2%AAAA, r2%ABAB, r2%BBBB, &
                              Wmnij%AAAA, Wmnij%ABAB, Wmnij%BBBB, fact, iCore, maxCore)
        else
            call c2_r2Wmnij_R(c2%ABAB, r2%ABAB, Wmnij%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2_r2Wmnij_SA error.'
            stop
        end if

    end subroutine c2_r2Wmnij_SA

    !> c(ab,ij) += 0.5*fact*r(ef,ij)*V(ab,ef)
    subroutine c2_r2abcd_SA(c2, r2, abcd, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r2, abcd
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abcd%spin_exist('ABAB')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abcd%spin_exist('AAAA')
            spin_exist = spin_exist .and. abcd%spin_exist('BBBB')

            call c2_r2abcd_U(c2%AAAA, c2%ABAB, c2%BBBB, r2%AAAA, r2%ABAB, r2%BBBB, &
                             abcd%AAAA, abcd%ABAB, abcd%BBBB, fact, iCore, maxCore)
        else
            call c2_r2abcd_R(c2%ABAB, r2%ABAB, abcd%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2_r2abcd_SA error.'
            stop
        end if

    end subroutine c2_r2abcd_SA

    !> c(ab,ij) -= fact*P(ab){r(a,m)*V(ij,mb)}
    subroutine c2_r1ijka_SA(c2, r1, ijka, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r1, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call c2_r1ijka_U(c2%AAAA, c2%ABAB, c2%BBBB, r1%AA, r1%BB, &
                             ijka%AAAA, ijka%ABAB, ijka%BBBB, ijka%BABA, fact, iCore, maxCore)
        else
            call c2_r1ijka_R(c2%ABAB, r1%AA, ijka%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2_r1ijka_SA error.'
            stop
        end if

    end subroutine c2_r1ijka_SA

    !> c(ab,ij) += fact*P(ab|ij){r(ae,im)*Wmbej(em,bj)}
    subroutine c2_r2Wmbej_SA(c2, r2, Wmbej, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r2, Wmbej
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Wmbej%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wmbej%spin_exist('AABB')
        spin_exist = spin_exist .and. Wmbej%spin_exist('ABAB')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBAA')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BABA')

            call c2_r2Wmbej_U(c2%AAAA, c2%ABAB, c2%BBBB, r2%AAAA, r2%ABAB, r2%BBBB, &
                              Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, Wmbej%BBBB, Wmbej%BBAA, Wmbej%BABA, &
                              fact, iCore, maxCore)
        else
            call c2_r2Wmbej_R(c2%ABAB, r2%AAAA, r2%ABAB, Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, &
                              fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2_r2Wmbej_SA error.'
            stop
        end if

    end subroutine c2_r2Wmbej_SA

    !> c(ab,ij) += fact*P(ij){r(e,i)*V(ab,ej)}
    subroutine c2_r1abci_SA(c2, r1, abci, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2
        type(SAopt), intent(in) :: r1, abci
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abci%spin_exist('ABAB')

        if (c2%UHF) then
            spin_exist = spin_exist .and. c2%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abci%spin_exist('AAAA')
            spin_exist = spin_exist .and. abci%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BABA')

            call c2_r1abci_U(c2%AAAA, c2%ABAB, c2%BBBB, r1%AA, r1%BB, &
                             abci%AAAA, abci%ABAB, abci%BBBB, abci%BABA, fact, iCore, maxCore)
        else
            call c2_r1abci_R(c2%ABAB, r1%AA, abci%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2_r1abci_SA error.'
            stop
        end if

    end subroutine c2_r1abci_SA

    !> F(a,e) -= fact*f(e,l)*r(a,l)
    subroutine Fae_r1fia_SA(Fae, r1, fia, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fae
        type(SAopt), intent(in) :: r1, fia
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fae%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. fia%spin_exist('AA')

        if (Fae%UHF) then
            spin_exist = spin_exist .and. Fae%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. fia%spin_exist('BB')

            call Fae_r1fia_U(Fae%AA, Fae%BB, r1%AA, r1%BB, fia%AA, fia%BB, &
                             fact, iCore, maxCore)
        else
            call Fae_r1fia_R(Fae%AA, r1%AA, fia%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fae_r1fia_SA error.'
            stop
        end if

    end subroutine Fae_r1fia_SA

    !> F(a,e) += fact*r(f,m)*V(ef,am)
    subroutine Fae_r1abci_SA(Fae, r1, abci, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fae
        type(SAopt), intent(in) :: r1, abci
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fae%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abci%spin_exist('AAAA')
        spin_exist = spin_exist .and. abci%spin_exist('ABAB')

        if (Fae%UHF) then
            spin_exist = spin_exist .and. Fae%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abci%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BABA')

            call Fae_r1abci_U(Fae%AA, Fae%BB, r1%AA, r1%BB, &
                              abci%AAAA, abci%ABAB, abci%BBBB, abci%BABA, &
                              fact, iCore, maxCore)
        else
            call Fae_r1abci_R(Fae%AA, r1%AA, abci%AAAA, abci%ABAB, &
                              fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fae_r1abci_SA error.'
            stop
        end if

    end subroutine Fae_r1abci_SA

    !> F(a,e) -= fact*0.5*r(af,mn)*V(ef,mn)
    subroutine Fae_r2abij_SA(Fae, r2, abij, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fae
        type(SAopt), intent(in) :: r2, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fae%spin_exist('AA')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abij%spin_exist('AAAA')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Fae%UHF) then
            spin_exist = spin_exist .and. Fae%spin_exist('BB')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Fae_r2abij_U(Fae%AA, Fae%BB, r2%AAAA, r2%ABAB, r2%BBBB, &
                              abij%AAAA, abij%ABAB, abij%BBBB, fact, iCore, maxCore)
        else
            call Fae_r2abij_R(Fae%AA, r2%AAAA, r2%ABAB, abij%AAAA, abij%ABAB, &
                              fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fae_r2abij_SA error.'
            stop
        end if

    end subroutine Fae_r2abij_SA

    !> F(m,i) += fact*r(d,i)*f(d,m)
    subroutine Fmi_r1fia_SA(Fmi, r1, fia, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fmi
        type(SAopt), intent(in) :: r1, fia
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fmi%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. fia%spin_exist('AA')

        if (Fmi%UHF) then
            spin_exist = spin_exist .and. Fmi%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. fia%spin_exist('BB')

            call Fmi_r1fia_U(Fmi%AA, Fmi%BB, r1%AA, r1%BB, fia%AA, fia%BB, &
                             fact, iCore, maxCore)
        else
            call Fmi_r1fia_R(Fmi%AA, r1%AA, fia%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fmi_r1fia_SA error.'
            stop
        end if

    end subroutine Fmi_r1fia_SA

    !> F(m,i) += fact*r(e,n)*V(mn,ie)
    subroutine Fmi_r1ijka_SA(Fmi, r1, ijka, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fmi
        type(SAopt), intent(in) :: r1, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fmi%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (Fmi%UHF) then
            spin_exist = spin_exist .and. Fmi%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call Fmi_r1ijka_U(Fmi%AA, Fmi%BB, r1%AA, r1%BB, &
                              ijka%AAAA, ijka%ABAB, ijka%BBBB, ijka%BABA, &
                              fact, iCore, maxCore)
        else
            call Fmi_r1ijka_R(Fmi%AA, r1%AA, ijka%AAAA, ijka%ABAB, &
                              fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fmi_r1ijka_SA error.'
            stop
        end if

    end subroutine Fmi_r1ijka_SA

    !> F(m,i) += fact*0.5*r(ef,in)*V(ef,mn)
    subroutine Fmi_r2abij_SA(Fmi, r2, abij, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fmi
        type(SAopt), intent(in) :: r2, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fmi%spin_exist('AA')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abij%spin_exist('AAAA')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Fmi%UHF) then
            spin_exist = spin_exist .and. Fmi%spin_exist('BB')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Fmi_r2abij_U(Fmi%AA, Fmi%BB, r2%AAAA, r2%ABAB, r2%BBBB, &
                              abij%AAAA, abij%ABAB, abij%BBBB, fact, iCore, maxCore)
        else
            call Fmi_r2abij_R(Fmi%AA, r2%AAAA, r2%ABAB, abij%AAAA, abij%ABAB, &
                              fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fmi_r2abij_SA error.'
            stop
        end if

    end subroutine Fmi_r2abij_SA

    !> F(e,m) += fact*r(f,n)*V(ef,mn)
    subroutine Fme_r1abij_SA(Fme, r1, abij, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Fme
        type(SAopt), intent(in) :: r1, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Fme%spin_exist('AA')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abij%spin_exist('AAAA')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Fme%UHF) then
            spin_exist = spin_exist .and. Fme%spin_exist('BB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Fme_r1abij_U(Fme%AA, Fme%BB, r1%AA, r1%BB, abij%AAAA, abij%ABAB, abij%BBBB, &
                              fact, iCore, maxCore)
        else
            call Fme_r1abij_R(Fme%AA, r1%AA, abij%AAAA, abij%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Fme_r1abij_SA error.'
            stop
        end if

    end subroutine Fme_r1abij_SA

    !> W(mn,ij) += fact*P(ij){r(e,j)*V(mn,ie)}
    subroutine Wmnij_r1ijka_SA(Wmnij, r1, ijka, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wmnij
        type(SAopt), intent(in) :: r1, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB

        character(len=3) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:3) = spin(1:3)
        else
            spin_case(1:3) = 'TTT'
        end if

        spin_exist = spin_exist .and. Wmnij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (Wmnij%UHF) then
            spin_exist = spin_exist .and. Wmnij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wmnij%spin_exist('BBBB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call Wmnij_r1ijka_U(Wmnij%AAAA, Wmnij%ABAB, Wmnij%BBBB, spin_case, &
                                r1%AA, r1%BB, ijka%AAAA, ijka%ABAB, ijka%BBBB, ijka%BABA, &
                                fact, iCore, maxCore)
        else
            call Wmnij_r1ijka_R(Wmnij%ABAB, r1%AA, ijka%ABAB, &
                                fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wmnij_r1ijka_SA error.'
            stop
        end if

    end subroutine Wmnij_r1ijka_SA

    !> W(mn,ij) += fact*0.5*r(ef,ij)*V(ef,mn)
    subroutine Wmnij_r2abij_SA(Wmnij, r2, abij, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wmnij
        type(SAopt), intent(in) :: r2, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB

        character(len=3) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:3) = spin(1:3)
        else
            spin_case(1:3) = 'TTT'
        end if

        spin_exist = spin_exist .and. Wmnij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Wmnij%UHF) then
            spin_exist = spin_exist .and. Wmnij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wmnij%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abij%spin_exist('AAAA')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Wmnij_r2abij_U(Wmnij%AAAA, Wmnij%ABAB, Wmnij%BBBB, &
                                r2%AAAA, r2%ABAB, r2%BBBB, abij%AAAA, abij%ABAB, abij%BBBB, &
                                fact, iCore, maxCore)
        else
            call Wmnij_r2abij_R(Wmnij%ABAB, r2%ABAB, abij%ABAB, &
                                fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wmnij_r2abij_SA error.'
            stop
        end if

    end subroutine Wmnij_r2abij_SA

    !> Wmbej(em,bj) += fact*r(f,j)*V(fe,bm)
    subroutine Wmbej_r1abci_SA(Wmbej, r1, abci, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wmbej
        type(SAopt), intent(in) :: r1, abci
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,AABB,ABAB,BBBB,BBAA,BABA

        character(len=6) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:6) = spin(1:6)
        else
            spin_case(1:6) = 'TTTTTT'
        end if

        spin_exist = spin_exist .and. Wmbej%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wmbej%spin_exist('AABB')
        spin_exist = spin_exist .and. Wmbej%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abci%spin_exist('AAAA')
        spin_exist = spin_exist .and. abci%spin_exist('ABAB')

        if (Wmbej%UHF) then
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBAA')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abci%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BABA')

            call Wmbej_r1abci_U(Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, Wmbej%BBBB, Wmbej%BBAA, Wmbej%BABA, spin_case, &
                                r1%AA, r1%BB, abci%AAAA, abci%ABAB, abci%BBBB, abci%BABA, &
                                fact, iCore, maxCore)
        else
            call Wmbej_r1abci_R(Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, spin_case(1:3), &
                                r1%AA, abci%AAAA, abci%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wmbej_r1abci_SA error.'
            stop
        end if

    end subroutine Wmbej_r1abci_SA

    !> Wmbej(em,bj) += fact*0.5*r(bf,jn)*V(fe,nm)
    subroutine Wmbej_r2abij_SA(Wmbej, r2, abij, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wmbej
        type(SAopt), intent(in) :: r2, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,AABB,ABAB,BBBB,BBAA,BABA

        character(len=6) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:6) = spin(1:6)
        else
            spin_case(1:6) = 'TTTTTT'
        end if

        spin_exist = spin_exist .and. Wmbej%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wmbej%spin_exist('AABB')
        spin_exist = spin_exist .and. Wmbej%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abij%spin_exist('AAAA')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Wmbej%UHF) then
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBAA')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Wmbej_r2abij_U(Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, Wmbej%BBBB, Wmbej%BBAA, Wmbej%BABA, spin_case, &
                                r2%AAAA, r2%ABAB, r2%BBBB, abij%AAAA, abij%ABAB, abij%BBBB, &
                                fact, iCore, maxCore)
        else
            call Wmbej_r2abij_R(Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, spin_case(1:3), &
                                r2%AAAA, r2%ABAB, abij%AAAA, abij%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wmbej_r2abij_SA error.'
            stop
        end if

    end subroutine Wmbej_r2abij_SA

    !> Wmbej(em,bj) -= fact*r(b,n)*V(nm,je)
    subroutine Wmbej_r1ijka_SA(Wmbej, r1, ijka, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wmbej
        type(SAopt), intent(in) :: r1, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,AABB,ABAB,BBBB,BBAA,BABA

        character(len=6) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:6) = spin(1:6)
        else
            spin_case(1:6) = 'TTTTTT'
        end if

        spin_exist = spin_exist .and. Wmbej%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wmbej%spin_exist('AABB')
        spin_exist = spin_exist .and. Wmbej%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (Wmbej%UHF) then
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BBAA')
            spin_exist = spin_exist .and. Wmbej%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call Wmbej_r1ijka_U(Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, Wmbej%BBBB, Wmbej%BBAA, Wmbej%BABA, spin_case, &
                                r1%AA, r1%BB, ijka%AAAA, ijka%ABAB, ijka%BABA, ijka%BBBB, &
                                fact, iCore, maxCore)
        else
            call Wmbej_r1ijka_R(Wmbej%AAAA, Wmbej%AABB, Wmbej%ABAB, spin_case(1:3), &
                                r1%AA, ijka%AAAA, ijka%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wmbej_r1ijka_SA error.'
            stop
        end if

    end subroutine Wmbej_r1ijka_SA

    !> W(ij,kb) += fact*0.5*V(dc,bk)*r(cd,ij), need V(dc,bk)=-V(cd,bk)
    subroutine Wkbij_r2abci_SA(Wkbij, r2, abci, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wkbij
        type(SAopt), intent(in) :: r2, abci
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wkbij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abci%spin_exist('ABAB')

        if (Wkbij%UHF) then
            spin_exist = spin_exist .and. Wkbij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('AAAA')
            spin_exist = spin_exist .and. abci%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BABA')

            call Wkbij_r2abci_U(Wkbij%AAAA, Wkbij%ABAB, Wkbij%BBBB, Wkbij%BABA, spin_case, &
                                r2%AAAA, r2%ABAB, r2%BBBB, abci%AAAA, abci%ABAB, abci%BBBB, abci%BABA, &
                                fact, iCore, maxCore)
        else
            call Wkbij_r2abci_R(Wkbij%ABAB, r2%ABAB, abci%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wkbij_r2abci_SA error.'
            stop
        end if

    end subroutine Wkbij_r2abci_SA

    !> W(ij,kb) -= fact*P(ij){Vkbcj(ck,bj)*r(c,i)}, need V(kb,cj)=V(jc,bk)
    subroutine Wkbij_r1iabj_SA(Wkbij, r1, iabj, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wkbij
        type(SAopt), intent(in) :: r1, iabj
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wkbij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. iabj%spin_exist('AABB')
        spin_exist = spin_exist .and. iabj%spin_exist('ABAB')

        if (Wkbij%UHF) then
            spin_exist = spin_exist .and. Wkbij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. iabj%spin_exist('AAAA')
            spin_exist = spin_exist .and. iabj%spin_exist('BBBB')
            spin_exist = spin_exist .and. iabj%spin_exist('BABA')

            call Wkbij_r1iabj_U(Wkbij%AAAA, Wkbij%ABAB, Wkbij%BBBB, Wkbij%BABA, spin_case, &
                                r1%AA, r1%BB, iabj%AAAA, iabj%AABB, iabj%ABAB, iabj%BBBB, iabj%BABA, &
                                fact, iCore, maxCore)
        else
            call Wkbij_r1iabj_R(Wkbij%ABAB, r1%AA, iabj%AABB, iabj%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wkbij_r1iabj_SA error.'
            stop
        end if

    end subroutine Wkbij_r1iabj_SA

    !> W(ij,kb) -= fact*F(c,k)*r(cb,ij)
    subroutine Wkbij_r2Fkc_SA(Wkbij, r2, Fkc, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wkbij
        type(SAopt), intent(in) :: r2, Fkc
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wkbij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fkc%spin_exist('AA')

        if (Wkbij%UHF) then
            spin_exist = spin_exist .and. Wkbij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fkc%spin_exist('BB')

            call Wkbij_r2Fkc_U(Wkbij%AAAA, Wkbij%ABAB, Wkbij%BBBB, Wkbij%BABA, spin_case, &
                               r2%AAAA, r2%ABAB, r2%BBBB, Fkc%AA, Fkc%BB, &
                               fact, iCore, maxCore)
        else
            call Wkbij_r2Fkc_R(Wkbij%ABAB, r2%ABAB, Fkc%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wkbij_r2Fkc_SA error.'
            stop
        end if

    end subroutine Wkbij_r2Fkc_SA

    !> W(ij,kb) -= fact*U(kl,ij)*r(b,l)
    subroutine Wkbij_r1Uklij_SA(Wkbij, r1, Uklij, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wkbij
        type(SAopt), intent(in) :: r1, Uklij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wkbij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. Uklij%spin_exist('ABAB')

        if (Wkbij%UHF) then
            spin_exist = spin_exist .and. Wkbij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. Uklij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Uklij%spin_exist('BBBB')

            call Wkbij_r1Uklij_U(Wkbij%AAAA, Wkbij%ABAB, Wkbij%BBBB, Wkbij%BABA, spin_case, &
                                 r1%AA, r1%BB, Uklij%AAAA, Uklij%ABAB, Uklij%BBBB, &
                                 fact, iCore, maxCore)
        else
            call Wkbij_r1Uklij_R(Wkbij%ABAB, r1%AA, Uklij%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wkbij_r1Uklij_SA error.'
            stop
        end if

    end subroutine Wkbij_r1Uklij_SA

    !> W(ij,kb) += fact*P(ij){V(kl,ic)*r(bc,jl)}
    subroutine Wkbij_r2ijka_SA(Wkbij, r2, ijka, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wkbij
        type(SAopt), intent(in) :: r2, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wkbij%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (Wkbij%UHF) then
            spin_exist = spin_exist .and. Wkbij%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wkbij%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call Wkbij_r2ijka_U(Wkbij%AAAA, Wkbij%ABAB, Wkbij%BBBB, Wkbij%BABA, spin_case, &
                                r2%AAAA, r2%ABAB, r2%BBBB, &
                                ijka%AAAA, ijka%ABAB, ijka%BBBB, ijka%BABA, &
                                fact, iCore, maxCore)
        else
            call Wkbij_r2ijka_R(Wkbij%ABAB, r2%AAAA, r2%ABAB, ijka%AAAA, ijka%ABAB, &
                                fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wkbij_r2ijka_SA error.'
            stop
        end if

    end subroutine Wkbij_r2ijka_SA

    !> W(lk,jd) += fact*r(c,j)*V(cd,lk)
    subroutine Wlkjd_r1abij_SA(Wlkjd, r1, abij, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wlkjd
        type(SAopt), intent(in) :: r1, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wlkjd%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wlkjd%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abij%spin_exist('AAAA')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Wlkjd%UHF) then
            spin_exist = spin_exist .and. Wlkjd%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wlkjd%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Wlkjd_r1abij_U(Wlkjd%AAAA, Wlkjd%ABAB, Wlkjd%BBBB, Wlkjd%BABA, spin_case, &
                                r1%AA, r1%BB, abij%AAAA, abij%ABAB, abij%BBBB, &
                                fact, iCore, maxCore)
        else
            call Wlkjd_r1abij_R(Wlkjd%AAAA, Wlkjd%ABAB, spin_case(1:2), &
                                r1%AA, abij%AAAA, abij%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wlkjd_r1abij_SA error.'
            stop
        end if

    end subroutine Wlkjd_r1abij_SA

    !> W(ab,dj) += fact*V(ab,dc)*r(c,j)
    subroutine Wabdj_r1abcd_SA(Wabdj, r1, abcd, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wabdj
        type(SAopt), intent(in) :: r1, abcd
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wabdj%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abcd%spin_exist('ABAB')

        if (Wabdj%UHF) then
            spin_exist = spin_exist .and. Wabdj%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BBBB')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abcd%spin_exist('AAAA')
            spin_exist = spin_exist .and. abcd%spin_exist('BBBB')

            call Wabdj_r1abcd_U(Wabdj%AAAA, Wabdj%ABAB, Wabdj%BBBB, Wabdj%BABA, spin_case, &
                                r1%AA, r1%BB, abcd%AAAA, abcd%ABAB, abcd%BBBB, &
                                fact, iCore, maxCore)
        else
            call Wabdj_r1abcd_R(Wabdj%ABAB, r1%AA, abcd%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wabdj_r1abcd_SA error.'
            stop
        end if

    end subroutine Wabdj_r1abcd_SA

    !> W(ab,dj) -= fact*P(ab){Ukbdj(dk,bj)*r(a,k)}
    subroutine Wabdj_r1Ukbdj_SA(Wabdj, r1, Ukbdj, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wabdj
        type(SAopt), intent(in) :: r1, Ukbdj
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wabdj%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. Ukbdj%spin_exist('AABB')
        spin_exist = spin_exist .and. Ukbdj%spin_exist('ABAB')

        if (Wabdj%UHF) then
            spin_exist = spin_exist .and. Wabdj%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. Ukbdj%spin_exist('AAAA')
            spin_exist = spin_exist .and. Ukbdj%spin_exist('BBBB')
            spin_exist = spin_exist .and. Ukbdj%spin_exist('BBAA')
            spin_exist = spin_exist .and. Ukbdj%spin_exist('BABA')

            call Wabdj_r1Ukbdj_U(Wabdj%AAAA, Wabdj%ABAB, Wabdj%BBBB, Wabdj%BABA, spin_case, &
                                 r1%AA, r1%BB, Ukbdj%AAAA, Ukbdj%AABB, Ukbdj%ABAB, &
                                 Ukbdj%BBBB, Ukbdj%BBAA, Ukbdj%BABA, &
                                 fact, iCore, maxCore)
        else
            call Wabdj_r1Ukbdj_R(Wabdj%ABAB, r1%AA, Ukbdj%AABB, Ukbdj%ABAB, &
                                 fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wabdj_r1Ukbdj_SA error.'
            stop
        end if

    end subroutine Wabdj_r1Ukbdj_SA

    !> W(ab,dj) += fact*0.5*V(lk,jd)*r(ab,kl)
    subroutine Wabdj_r2ijka_SA(Wabdj, r2, ijka, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wabdj
        type(SAopt), intent(in) :: r2, ijka
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wabdj%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. ijka%spin_exist('ABAB')

        if (Wabdj%UHF) then
            spin_exist = spin_exist .and. Wabdj%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('AAAA')
            spin_exist = spin_exist .and. ijka%spin_exist('BBBB')
            spin_exist = spin_exist .and. ijka%spin_exist('BABA')

            call Wabdj_r2ijka_U(Wabdj%AAAA, Wabdj%ABAB, Wabdj%BBBB, Wabdj%BABA, spin_case, &
                                r2%AAAA, r2%ABAB, r2%BBBB, &
                                ijka%AAAA, ijka%ABAB, ijka%BBBB, ijka%BABA, &
                                fact, iCore, maxCore)
        else
            call Wabdj_r2ijka_R(Wabdj%ABAB, r2%ABAB, ijka%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wabdj_r2ijka_SA error.'
            stop
        end if

    end subroutine Wabdj_r2ijka_SA

    !> W(ab,dj) -= fact*F(d,l)*r(ab,lj)
    subroutine Wabdj_r2fld_SA(Wabdj, r2, Fld, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wabdj
        type(SAopt), intent(in) :: r2, Fld
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wabdj%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fld%spin_exist('AA')

        if (Wabdj%UHF) then
            spin_exist = spin_exist .and. Wabdj%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fld%spin_exist('BB')

            call Wabdj_r2Fld_U(Wabdj%AAAA, Wabdj%ABAB, Wabdj%BBBB, Wabdj%BABA, spin_case, &
                               r2%AAAA, r2%ABAB, r2%BBBB, Fld%AA, Fld%BB, &
                               fact, iCore, maxCore)
        else
            call Wabdj_r2Fld_R(Wabdj%ABAB, r2%ABAB, Fld%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wabdj_r2fld_SA error.'
            stop
        end if

    end subroutine Wabdj_r2fld_SA

    !> W(ab,dj) += fact*P(ab){V(dc,ak)*r(bc,jk)}
    subroutine Wabdj_r2abci_SA(Wabdj, r2, abci, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wabdj
        type(SAopt), intent(in) :: r2, abci
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wabdj%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('AAAA')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')
        spin_exist = spin_exist .and. abci%spin_exist('AAAA')
        spin_exist = spin_exist .and. abci%spin_exist('ABAB')

        if (Wabdj%UHF) then
            spin_exist = spin_exist .and. Wabdj%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wabdj%spin_exist('BABA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BBBB')
            spin_exist = spin_exist .and. abci%spin_exist('BABA')

            call Wabdj_r2abci_U(Wabdj%AAAA, Wabdj%ABAB, Wabdj%BBBB, Wabdj%BABA, spin_case, &
                                r2%AAAA, r2%ABAB, r2%BBBB, &
                                abci%AAAA, abci%ABAB, abci%BBBB, abci%BABA, &
                                fact, iCore, maxCore)
        else
            call Wabdj_r2abci_R(Wabdj%ABAB, r2%AAAA, r2%ABAB, abci%AAAA, abci%ABAB, &
                                fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wabdj_r2abci_SA error.'
            stop
        end if

    end subroutine Wabdj_r2abci_SA

    !> W(cd,ak) -= fact*V(cd,lk)*r(a,l)
    subroutine Wakcd_r1abij_SA(Wakcd, r1, abij, fact, iCore, maxCore, spin)
        implicit none

        type(SAopt), intent(inout) :: Wakcd
        type(SAopt), intent(in) :: r1, abij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore
        character(len=*), intent(in), optional :: spin !> AAAA,ABAB,BBBB,BABA

        character(len=4) :: spin_case
        logical :: spin_exist = T

        if (present(spin)) then
            spin_case(1:4) = spin(1:4)
        else
            spin_case(1:4) = 'TTTT'
        end if

        spin_exist = spin_exist .and. Wakcd%spin_exist('ABAB')
        spin_exist = spin_exist .and. r1%spin_exist('AA')
        spin_exist = spin_exist .and. abij%spin_exist('ABAB')

        if (Wakcd%UHF) then
            spin_exist = spin_exist .and. Wakcd%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wakcd%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wakcd%spin_exist('BABA')
            spin_exist = spin_exist .and. r1%spin_exist('BB')
            spin_exist = spin_exist .and. abij%spin_exist('AAAA')
            spin_exist = spin_exist .and. abij%spin_exist('BBBB')

            call Wakcd_r1abij_U(Wakcd%AAAA, Wakcd%ABAB, Wakcd%BBBB, Wakcd%BABA, spin_case, &
                                r1%AA, r1%BB, abij%AAAA, abij%ABAB, abij%BBBB, &
                                fact, iCore, maxCore)
        else
            call Wakcd_r1abij_R(Wakcd%ABAB, r1%AA, abij%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wakcd_r1abij_SA error.'
            stop
        end if

    end subroutine Wakcd_r1abij_SA

    !> c(a,i) += fact*f(i,j)*l(a,j)
    subroutine c1d_l1fij_SA(c1d, l1, fij, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: l1, fij
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1d%spin_exist('AA')
        spin_exist = spin_exist .and. l1%spin_exist('AA')
        spin_exist = spin_exist .and. fij%spin_exist('AA')

        if (c1d%UHF) then
            spin_exist = spin_exist .and. c1d%spin_exist('BB')
            spin_exist = spin_exist .and. l1%spin_exist('BB')
            spin_exist = spin_exist .and. fij%spin_exist('BB')

            call c1d_l1fij_U(c1d%AA, c1d%BB, l1%AA, l1%BB, &
                             fij%AA, fij%BB, fact, iCore, maxCore)
        else
            call c1d_l1fij_R(c1d%AA, l1%AA, fij%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1d_l1fij_SA error.'
            stop
        end if

    end subroutine c1d_l1fij_SA

    !> c(a,i) -= fact*f(b,a)*l(b,i)
    subroutine c1d_l1fab_SA(c1d, l1, fab, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: l1, fab
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1d%spin_exist('AA')
        spin_exist = spin_exist .and. l1%spin_exist('AA')
        spin_exist = spin_exist .and. fab%spin_exist('AA')

        if (c1d%UHF) then
            spin_exist = spin_exist .and. c1d%spin_exist('BB')
            spin_exist = spin_exist .and. l1%spin_exist('BB')
            spin_exist = spin_exist .and. fab%spin_exist('BB')

            call c1d_l1fab_U(c1d%AA, c1d%BB, l1%AA, l1%BB, &
                             fab%AA, fab%BB, fact, iCore, maxCore)
        else
            call c1d_l1fab_R(c1d%AA, l1%AA, fab%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1d_l1fab_SA error.'
            stop
        end if

    end subroutine c1d_l1fab_SA

    !> c(a,i) += fact*Wibaj(ai,bj)*l(b,j)
    subroutine c1d_l1Wibaj_SA(c1d, l1, Wibaj, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: l1, Wibaj
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1d%spin_exist('AA')
        spin_exist = spin_exist .and. l1%spin_exist('AA')
        spin_exist = spin_exist .and. Wibaj%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wibaj%spin_exist('AABB')

        if (c1d%UHF) then
            spin_exist = spin_exist .and. c1d%spin_exist('BB')
            spin_exist = spin_exist .and. l1%spin_exist('BB')
            spin_exist = spin_exist .and. Wibaj%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wibaj%spin_exist('BBAA')

            call c1d_l1Wibaj_U(c1d%AA, c1d%BB, l1%AA, l1%BB, &
                               Wibaj%AAAA, Wibaj%AABB, Wibaj%BBBB, Wibaj%BBAA, &
                               fact, iCore, maxCore)
        else
            call c1d_l1Wibaj_R(c1d%AA, l1%AA, Wibaj%AAAA, Wibaj%AABB, &
                               fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1d_l1Wibaj_SA error.'
            stop
        end if

    end subroutine c1d_l1Wibaj_SA

    !> c(a,i) -= 0.5*fact*l(ab,kj)*W(kj,ib)
    !> please use c1_r2ijka_SA
    subroutine c1d_l2Wibkj_SA(c1d, l2, Wibkj, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: l2, Wibkj
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        call c1_r2ijka_SA(c1d, l2, Wibkj, fact, iCore, maxCore)

    end subroutine c1d_l2Wibkj_SA

    !> c(a,i) += fact*0.5*W(cb,aj)*l(cb,ij)
    !> please use c1_r2abci_SA
    subroutine c1d_l2Wcbaj_SA(c1d, l2, Wcbaj, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: l2, Wcbaj
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        call c1_r2abci_SA(c1d, l2, Wcbaj, fact, iCore, maxCore)

    end subroutine c1d_l2Wcbaj_SA

    !> c(a,i) -= fact*W(ca,di)*F(c,d)
    subroutine c1d_FcdWdica_SA(c1d, Fcd, Wdica, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: Fcd, Wdica
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1d%spin_exist('AA')
        spin_exist = spin_exist .and. Fcd%spin_exist('AA')
        spin_exist = spin_exist .and. Wdica%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wdica%spin_exist('ABAB')

        if (c1d%UHF) then
            spin_exist = spin_exist .and. c1d%spin_exist('BB')
            spin_exist = spin_exist .and. Fcd%spin_exist('BB')
            spin_exist = spin_exist .and. Wdica%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wdica%spin_exist('BABA')

            call c1d_FcdWdica_U(c1d%AA, c1d%BB, Fcd%AA, Fcd%BB, &
                                Wdica%AAAA, Wdica%ABAB, Wdica%BBBB, Wdica%BABA, &
                                fact, iCore, maxCore)
        else
            call c1d_FcdWdica_R(c1d%AA, Fcd%AA, Wdica%AAAA, Wdica%ABAB, &
                                fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1d_FcdWdica_SA error.'
            stop
        end if

    end subroutine c1d_FcdWdica_SA

    !> c(a,i) -= fact*W(li,ka)*F(k,l)
    subroutine c1d_FklWlika_SA(c1d, Fkl, Wlika, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c1d
        type(SAopt), intent(in) :: Fkl, Wlika
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c1d%spin_exist('AA')
        spin_exist = spin_exist .and. Fkl%spin_exist('AA')
        spin_exist = spin_exist .and. Wlika%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wlika%spin_exist('ABAB')

        if (c1d%UHF) then
            spin_exist = spin_exist .and. c1d%spin_exist('BB')
            spin_exist = spin_exist .and. Fkl%spin_exist('BB')
            spin_exist = spin_exist .and. Wlika%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wlika%spin_exist('BABA')

            call c1d_FklWlika_U(c1d%AA, c1d%BB, Fkl%AA, Fkl%BB, &
                                Wlika%AAAA, Wlika%ABAB, Wlika%BBBB, Wlika%BABA, &
                                fact, iCore, maxCore)
        else
            call c1d_FklWlika_R(c1d%AA, Fkl%AA, Wlika%AAAA, Wlika%ABAB, &
                                fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c1d_FklWlika_SA error.'
            stop
        end if

    end subroutine c1d_FklWlika_SA

    !> c(ab,ij) -= fact*P(ij){F(j,k)*l(ab,ik)}
    subroutine c2d_l2Fjk_SA(c2d, l2, Fjk, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: l2, Fjk
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2d%spin_exist('ABAB')
        spin_exist = spin_exist .and. l2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fjk%spin_exist('AA')

        if (c2d%UHF) then
            spin_exist = spin_exist .and. c2d%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2d%spin_exist('BBBB')
            spin_exist = spin_exist .and. l2%spin_exist('AAAA')
            spin_exist = spin_exist .and. l2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fjk%spin_exist('BB')

            call c2d_l2Fjk_U(c2d%AAAA, c2d%ABAB, c2d%BBBB, l2%AAAA, l2%ABAB, l2%BBBB, &
                             Fjk%AA, Fjk%BB, fact, iCore, maxCore)
        else
            call c2d_l2Fjk_R(c2d%ABAB, l2%ABAB, Fjk%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2d_l2Fjk_SA error.'
            stop
        end if

    end subroutine c2d_l2Fjk_SA

    !> c(ab,ij) += fact*P(ab){F(c,b)*l(ac,ij)}
    subroutine c2d_l2Fcb_SA(c2d, l2, Fcb, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: l2, Fcb
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2d%spin_exist('ABAB')
        spin_exist = spin_exist .and. l2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Fcb%spin_exist('AA')

        if (c2d%UHF) then
            spin_exist = spin_exist .and. c2d%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2d%spin_exist('BBBB')
            spin_exist = spin_exist .and. l2%spin_exist('AAAA')
            spin_exist = spin_exist .and. l2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Fcb%spin_exist('BB')

            call c2d_l2Fcb_U(c2d%AAAA, c2d%ABAB, c2d%BBBB, l2%AAAA, l2%ABAB, l2%BBBB, &
                             Fcb%AA, Fcb%BB, fact, iCore, maxCore)
        else
            call c2d_l2Fcb_R(c2d%ABAB, l2%ABAB, Fcb%AA, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2d_l2Fcb_SA error.'
            stop
        end if

    end subroutine c2d_l2Fcb_SA

    !> c(ab,ij) -= fact*P(ab){W(ij,kb)*l(a,k)}
    !> please use c2_r1ijka_SA
    subroutine c2d_l1Wijkb_SA(c2d, l1, Wijkb, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: l1, Wijkb
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        call c2_r1ijka_SA(c2d, l1, Wijkb, fact, iCore, maxCore)

    end subroutine c2d_l1Wijkb_SA

    !> c(ab,ij) += fact*P(ij){W(ab,cj)*l(c,i)}
    !> pleasu use c2_r1abci_SA
    subroutine c2d_l1Wciab_SA(c2d, l1, Wciab, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: l1, Wciab
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        call c2_r1abci_SA(c2d, l1, Wciab, fact, iCore, maxCore)

    end subroutine c2d_l1Wciab_SA

    !> c(ab,ij) += fact*P(ab|ij){Wjcbk(bj,ck)*l(ac,ik)}
    subroutine c2d_l2Wjcbk_SA(c2d, l2, Wjcbk, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: l2, Wjcbk
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2d%spin_exist('ABAB')
        spin_exist = spin_exist .and. l2%spin_exist('AAAA')
        spin_exist = spin_exist .and. l2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Wjcbk%spin_exist('AAAA')
        spin_exist = spin_exist .and. Wjcbk%spin_exist('AABB')
        spin_exist = spin_exist .and. Wjcbk%spin_exist('ABAB')

        if (c2d%UHF) then
            spin_exist = spin_exist .and. c2d%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2d%spin_exist('BBBB')
            spin_exist = spin_exist .and. l2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wjcbk%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wjcbk%spin_exist('BBAA')
            spin_exist = spin_exist .and. Wjcbk%spin_exist('BABA')

            call c2d_l2Wjcbk_U(c2d%AAAA, c2d%ABAB, c2d%BBBB, &
                               l2%AAAA, l2%ABAB, l2%BBBB, &
                               Wjcbk%AAAA, Wjcbk%AABB, Wjcbk%ABAB, &
                               Wjcbk%BBBB, Wjcbk%BBAA, Wjcbk%BABA, &
                               fact, iCore, maxCore)
        else
            call c2d_l2Wjcbk_R(c2d%ABAB, l2%AAAA, l2%ABAB, &
                               Wjcbk%AAAA, Wjcbk%AABB, Wjcbk%ABAB, &
                               fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2d_l2Wjcbk_SA error.'
            stop
        end if

    end subroutine c2d_l2Wjcbk_SA

    !> c(ab,ij) += fact*0.5*W(ij,kl)*l(ab,kl)
    subroutine c2d_l2Wijkl_SA(c2d, l2, Wijkl, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: l2, Wijkl
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2d%spin_exist('ABAB')
        spin_exist = spin_exist .and. l2%spin_exist('ABAB')
        spin_exist = spin_exist .and. Wijkl%spin_exist('ABAB')

        if (c2d%UHF) then
            spin_exist = spin_exist .and. c2d%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2d%spin_exist('BBBB')
            spin_exist = spin_exist .and. l2%spin_exist('AAAA')
            spin_exist = spin_exist .and. l2%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wijkl%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wijkl%spin_exist('BBBB')

            call c2d_l2Wijkl_U(c2d%AAAA, c2d%ABAB, c2d%BBBB, l2%AAAA, l2%ABAB, l2%BBBB, &
                               Wijkl%AAAA, Wijkl%ABAB, Wijkl%BBBB, fact, iCore, maxCore)
        else
            call c2d_l2Wijkl_R(c2d%ABAB, l2%ABAB, Wijkl%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2d_l2Wijkl_SA error.'
            stop
        end if

    end subroutine c2d_l2Wijkl_SA

    !> c(ab,ij) -= fact*W(ij,kd)*U(ba,dk)
    subroutine c2d_WijkdUdkba_SA(c2d, Wijkd, Udkba, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: c2d
        type(SAopt), intent(in) :: Wijkd, Udkba
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. c2d%spin_exist('ABAB')
        spin_exist = spin_exist .and. Wijkd%spin_exist('ABAB')
        spin_exist = spin_exist .and. Udkba%spin_exist('ABAB')

        if (c2d%UHF) then
            spin_exist = spin_exist .and. c2d%spin_exist('AAAA')
            spin_exist = spin_exist .and. c2d%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wijkd%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wijkd%spin_exist('BBBB')
            spin_exist = spin_exist .and. Wijkd%spin_exist('BABA')
            spin_exist = spin_exist .and. Udkba%spin_exist('AAAA')
            spin_exist = spin_exist .and. Udkba%spin_exist('BBBB')
            spin_exist = spin_exist .and. Udkba%spin_exist('BABA')

            call c2d_WijkdUdkba_U(c2d%AAAA, c2d%ABAB, c2d%BBBB, &
                                  Wijkd%AAAA, Wijkd%ABAB, Wijkd%BBBB, Wijkd%BABA, &
                                  Udkba%AAAA, Udkba%ABAB, Udkba%BBBB, Udkba%BABA, &
                                  fact, iCore, maxCore)
        else
            call c2d_WijkdUdkba_R(c2d%ABAB, Wijkd%ABAB, Udkba%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in c2d_WijkdUdkba_SA error.'
            stop
        end if

    end subroutine c2d_WijkdUdkba_SA

    !> W(ab,cd) += fact*1/2*l(ab,ij)*r(cd,ij)
    subroutine Wabcd_l2r2_SA(Wabcd, l2, r2, fact, iCore, maxCore)
        implicit none

        type(SAopt), intent(inout) :: Wabcd
        type(SAopt), intent(in) :: l2, r2
        real*8, intent(in) :: fact
        real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
        integer*8, intent(in) :: maxCore

        logical :: spin_exist = T

        spin_exist = spin_exist .and. Wabcd%spin_exist('ABAB')
        spin_exist = spin_exist .and. l2%spin_exist('ABAB')
        spin_exist = spin_exist .and. r2%spin_exist('ABAB')

        if (Wabcd%UHF) then
            spin_exist = spin_exist .and. Wabcd%spin_exist('AAAA')
            spin_exist = spin_exist .and. Wabcd%spin_exist('BBBB')
            spin_exist = spin_exist .and. l2%spin_exist('AAAA')
            spin_exist = spin_exist .and. l2%spin_exist('BBBB')
            spin_exist = spin_exist .and. r2%spin_exist('AAAA')
            spin_exist = spin_exist .and. r2%spin_exist('BBBB')

            call Wabcd_l2r2_U(Wabcd%AAAA, Wabcd%ABAB, Wabcd%BBBB, &
                              l2%AAAA, l2%ABAB, l2%BBBB, &
                              r2%AAAA, r2%ABAB, r2%BBBB, &
                              fact, iCore, maxCore)
        else
            call Wabcd_l2r2_R(Wabcd%ABAB, l2%ABAB, r2%ABAB, fact, iCore, maxCore)
        end if

        if (.not. spin_exist) then
            write (6, *) 'ERROR: spin in Wabcd_l2r2_SA error.'
            stop
        end if

    end subroutine Wabcd_l2r2_SA

end module mod_mult
