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
! MODULE hard_chains_mod
!  This module contains the hard-chain subroutines for the PC-SAFT equation
!  of state.
!
!   - pure real function a_tilde_hc(saft_para, rho, Temp, x)
!     ~ <type(saft_parameter)>    :: saft_para      ! PC-SAFT parameters
!     ~ <real, dimension(N_comp)> :: rho            ! densities
!     ~ <real, dimension(N_comp)> :: x              ! mole fractions
!
!
!   - subroutine ...
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module hard_chains_mod
  use kinds_mod, only: dp                           ! loading 'real64' precision specifier (compiler independent)
  use saft_parameter_mod, only: saft_parameter      ! loading derived type 'saft_parameter'
  use hard_spheres_mod, only: a_tilde_hs, g_hs_ij   ! loading hard-sphere Helmholtz energy density and radial distribution function
  implicit none
  private
  public :: a_tilde_hc

contains


  ! ############################################################################
  ! Function for the calculation of the hard-chain contribution to the reduced
  !  Helmholtz energy 'a_tilde_hc = A_hc/NkT' according to equation (A.) in [1]
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}

  ! pure real function (in<PC-SAFT parameter>, in<densites>, in<mole fractions>)
  !  in  < saft_para    :: saft_parameter {d(N_Comp)} >
  !  in  < rho          :: real(N_comp) >
  !  in  < x            :: real(N_comp) >
  !  out < a_tilde_hc   :: real         >
  ! ----------------------------------------------------------------------------
  pure real(dp) function a_tilde_hc(saft_para, rho, x)
    ! Input variables
    type(saft_parameter), intent(in)    :: saft_para    ! parameter container type
    real(dp), intent(in)                :: rho          ! densities
    real(dp), dimension(:), intent(in)  :: x            ! molar fractions

    ! Local variables
    integer   :: i
    real(dp)  :: m_dash    ! mean segment number
    real(dp)  :: summand

    ! Calculate equation (A.4) in [1]
    summand = 0.0_dp
    do i = 1, saft_para%N_comp
      summand = summand + (x(i)*saft_para%m(i) - 1._dp) &
              & * log(g_hs_ij(saft_para, rho, x, i, i))
    end do
    ! Mean segment number
    m_dash = sum(saft_para%m*x)                                     ! equation (A.5) in [1] 

    ! Calculate Helmholtz energy density for hard-chains
    a_tilde_hc = m_dash*a_tilde_hs(saft_para, rho, x) - summand     ! equation (A.4) in [1]
  end function a_tilde_hc
  ! ############################################################################



end module hard_chains_mod
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
