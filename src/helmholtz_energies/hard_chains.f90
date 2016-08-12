!  ______                                              ______  _______  _______  _______
! (______)                                            / _____)(_______)(_______)(_______)
!  _     _   ____  _____   ____   ___   ____   _____ ( (____   _______  _____       _
! | |   | | / ___)(____ | / _  | / _ \ |  _ \ (_____) \____ \ |  ___  ||  ___)     | |
! | |__/ / | |    / ___ |( (_| || |_| || | | |        _____) )| |   | || |         | |
! |_____/  |_|    \_____| \___ | \___/ |_| |_|       (______/ |_|   |_||_|         |_|
!                       (_____|
!
!  Written by: Gernot Bauer, Rolf Stierle
!  Last update: 7. August 2016
!
!  To-Do's:
!   - Comments
!   -


module hard_chains_mod
    use kinds_mod, only: dp
    use saft_parameter_mod, only: saft_parameter
    use hard_spheres_mod, only: a_tilde_hs, g_hs_ij
    implicit none
    private
    public :: a_tilde_hc

contains

    pure real(dp) function a_tilde_hc(saft_para, rho, x)
        type(saft_parameter), intent(in) :: saft_para ! parameter container type
        real(dp), intent(in) :: rho ! density
        real(dp), dimension(:), intent(in) :: x ! molefraction

        integer :: i
        real(dp) :: m_dash ! mean segment number
        real(dp) :: summand

        summand = 0.0_dp
        do i = 1, saft_para%N_comp
            summand = summand + (x(i)*saft_para%m(i) - 1._dp) &
                              * log(g_hs_ij(saft_para, rho, x, i, i))
        end do
        m_dash = sum(saft_para%m*x)
        a_tilde_hc = m_dash*a_tilde_hs(saft_para, rho, x) - summand
    end function a_tilde_hc
end module hard_chains_mod
