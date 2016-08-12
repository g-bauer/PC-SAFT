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

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! MODULE hard_spheres_mod
!  This module contains all dispersion contributions for the PC-SAFT equation
!  of state:
!
!   - pure real function a_tilde_hs(saft_para, rho, Temp, x)
!     ~ <type(saft_parameter)>    :: saft_para      ! PC-SAFT parameters
!     ~ <real, dimension(N_comp)> :: rho            ! densities
!     ~ <real, dimension(N_comp)> :: x              ! mole fractions
!
!   - subroutine ....
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module hard_spheres_mod
  use kinds_mod, only: dp                                 ! loading 'real64' precision specifier (compiler independent)
  use parameter_mod, only: PI                             ! loading pi=3.14159
  use saft_parameter_mod, only: saft_parameter            ! loading derived type 'saft_parameter'
  implicit none
  private                                                 ! declaring everythin in this module 'private'
  public :: a_tilde_hs                                    ! specifically declaring 'public'
  public :: g_hs_ij                                       ! specifically declaring 'public'

contains

  ! ############################################################################
  ! Function for the calculation of the hard-sphere contribution to the reduced
  !  Helmholtz energy 'a_tilde_hs = A_hs/NkT' according to equation (A.6) in [1]
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}

  ! pure real function (in<PC-SAFT parameter>, in<densites>, in<mole fractions>)
  !  in  < saft_para    :: saft_parameter {d(N_Comp)} >
  !  in  < rho          :: real(N_comp) >
  !  in  < x            :: real(N_comp) >
  !  out < a_tilde_disp :: real         >
  ! ----------------------------------------------------------------------------
  pure real(dp) function a_tilde_hs(saft_para, rho, x)
    ! Input parameters
    type(saft_parameter), intent(in)   :: saft_para   ! parameter container type
    real(dp), intent(in)               :: rho         ! density
    real(dp), dimension(:), intent(in) :: x           ! molefraction

    ! Local parameters
    real(dp), dimension(0:3) :: zeta ! weighted densities
    integer :: n


    ! Calculation of the weighted densities
    do n = 0, 3
      zeta(n) = PI/6._dp * rho * sum(x * saft_para%m * saft_para%d**n)                 ! equation (9) or (A.8) in [1]
    end do

    ! Reduced Helmholtz energy density contribution of the hard-sphere systems
    a_tilde_hs = 1/zeta(0) * ( 3 * zeta(1) * zeta(2) / (1-zeta(3)) &                   ! equation (A.6) in [1]
               & + zeta(2)**3 / (zeta(3) * (1-zeta(3))**2) &
               & + (zeta(2)**3 / zeta(3)**2 - zeta(0)) * log(1-zeta(3)) )
  end function a_tilde_hs
  ! ############################################################################


  ! ############################################################################
  ! Function for the calculation of the hard-sphere radial distribution function
  !  'g_hs(r)' according to equations (8) or (A.7) in [1].
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}
  !
  ! pure real function (in<PC-SAFT parameter>, in<densities>, in<molar fractions>, in<component index>, in<component index>)
  !  in  < saft_para :: saft_parameter {m(N_comp), d(N_comp)}>
  !  in  < rho       :: real(N_comp) >
  !  in  < x         :: real(N_comp) >
  !  in  < i         :: integer      >
  !  in  < j         :: integer      >
  ! ----------------------------------------------------------------------------------
  pure real(dp) function g_hs_ij(saft_para, rho, x, i, j)
    ! Input variables
    type(saft_parameter), intent(in)    :: saft_para    ! parameter container type
    real(dp), intent(in)                :: rho          ! density
    real(dp), dimension(:), intent(in)  :: x            ! molefraction
    integer, intent(in)                 :: i, j         ! component identifier

    ! Local variables
    real(dp), dimension(0:3)  :: zeta   ! weighted densities
    integer                   :: n      ! iteration variable


    ! Weighted densities
    do n = 0, 3
      zeta(n) = PI/6._dp * rho * sum(x * saft_para%m * saft_para%d**n)    ! equation (9) or (A.8) in [1]
    end do

    ! Radial distribution function 'g_ij(r)' for the hard-sphere fluid;
    !  equation (8) or (A.7) in [1]
    g_hs_ij = 1._dp / (1._dp-zeta(3))             &
            & + (saft_para%d(i) * saft_para%d(j)) / (saft_para%d(i) + saft_para%d(j)) &
            & * 3._dp*zeta(2) / (1._dp-zeta(3))**2 &
            & + ((saft_para%d(i) * saft_para%d(j)) / (saft_para%d(i) + saft_para%d(j)))**2 &
            & * zeta(2)**2 / (1._dp-zeta(3))**3
  end function g_hs_ij
  ! ############################################################################




end module hard_spheres_mod
