!  ______                                    ______ _______ _______ _______
! (______)                                  / _____|_______|_______|_______)
!  _     _ ____ _____  ____  ___  ____ ____( (____  _______ _____      _
! | |   | / ___|____ |/ _  |/ _ \|  _ (_____)____ \|  ___  |  ___)    | |
! | |__/ / |   / ___ ( (_| | |_| | | | |    _____) ) |   | | |        | |
! |_____/|_|   \_____|\___ |\___/|_| |_|   (______/|_|   |_|_|        |_|
!                    (_____|
!
!  Written by: Gernot Bauer, Rolf Stierle
!  Last update: 7. August 2016
!
!  To-Do's:
!   - Comments
!   -

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! MODULE hard_spheres_mod
!  This module contains all dispersion contributions for the PC-SAFT equation
!  of state:
!   - pure real function a_tilde_hs(saft_para, rho, Temp, x)
!     ~ <type(saft_parameter)>    :: saft_para      ! PC-SAFT parameters
!     ~ <real, dimension(N_comp)> :: rho            ! densities
!     ~ <real>                    :: Temp           ! temperature
!     ~ <real, dimension(N_comp)> :: x              ! mole fractions
!
!   - subroutine ....
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module hard_spheres_mod
    use kinds_mod, only: dp
    use parameter_mod, only: PI
    use saft_parameter_mod, only: saft_parameter
    implicit none
    private
    public :: a_tilde_hs
    public :: g_hs_ij

contains

    ! Contribution of the hard-sphere fluid to the Helmholtz energy a_tilde_hs = A_hs/NkT
    ! According to Equation (A.6) in [1]
    pure real(dp) function a_tilde_hs(saft_para, rho, x)
      type(saft_parameter), intent(in) :: saft_para ! parameter container type
      real(dp), intent(in) :: rho ! density
      real(dp), dimension(:), intent(in) :: x ! molefraction

      real(dp), dimension(0:3) :: zeta ! weighted densities
      integer :: n

      do n = 0, 3
        zeta(n) = PI/6._dp * rho * sum(x * saft_para%m * saft_para%d**n)
      end do

      a_tilde_hs = 1/zeta(0) * ( 3 * zeta(1) * zeta(2) / (1-zeta(3)) &
                 & + zeta(2)**3 / (zeta(3) * (1-zeta(3))**2) &
                 & + (zeta(2)**3 / zeta(3)**2 - zeta(0)) * log(1-zeta(3)) )

    end function a_tilde_hs



    !
    pure real(dp) function g_hs_ij(saft_para, rho, x, i, j)
      type(saft_parameter), intent(in) :: saft_para ! parameter container type
      real(dp), intent(in) :: rho ! density
      real(dp), dimension(:), intent(in) :: x ! molefraction
      integer, intent(in) :: i, j

      real(dp), dimension(0:3) :: zeta ! weighted densities
      integer :: n

      do n = 0, 3
        zeta(n) = PI/6._dp * rho * sum(x * saft_para%m * saft_para%d**n)
      end do

      g_hs_ij = 1._dp / (1._dp-zeta(3)) + (saft_para%d(i) * saft_para%d(j)) / (saft_para%d(i) + saft_para%d(j)) &
           & * 3._dp*zeta(2) / (1._dp-zeta(3))**2 &
           & + ((saft_para%d(i) * saft_para%d(j)) / (saft_para%d(i) + saft_para%d(j)))**2 &
           & * zeta(2)**2 / (1._dp-zeta(3))**3
    end function g_hs_ij




end module hard_spheres_mod
