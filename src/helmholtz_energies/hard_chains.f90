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
  use hard_spheres_mod, only: a_tilde_hs, g_hs_ij, Z_hs, rho_dg_hs_ij_drho
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





  ! ############################################################################
  ! Function for the calculation of the hard-chain compressibility 'Z_hc'
  !  according to equations (A.25) in [1].
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}
  !
  ! pure real function (in<PC-SAFT parameter>, in<densities>, in<molar fractions>)
  !  in  < saft_para :: saft_parameter {m(N_comp), d(N_comp)}>
  !  in  < rho       :: real(N_comp) >
  !  in  < x         :: real(N_comp) >
  ! ----------------------------------------------------------------------------
  pure real(dp) function Z_hc(saft_para, rho, x)
    ! Input variables
    type(saft_parameter)                  :: saft_para    ! parameter container type
    real(dp), dimension(:,), intent(in)   :: rho          ! density
    real(dp), intent(in)                  :: x            ! molefraction

    ! Local variables
    real(dp)                  :: m_dash     ! mean segment number
    real(dp)                  :: summand
    integer                   :: i          ! iteration variable

    ! Mean segment number
    m_dash = sum(saft_para%m*x)                                     ! equation (A.5) in [1]

    !
    summand = 0.0
    do i = 1, saft_para%N_comp
      summand = summand + x(i)*(m_dash - 1._dp)/g_hs_ij(saft_para,rho,x,i,i) &  ! equation (A.25) in [1]
              & * rho_dg_hs_ij_drho(saft_para,rho,x,i,i)
    end do


    ! Calculation of compressibility 'Z_hs'
    Z_hc = m_dash*Z_hs(saft_para,rho_x) - summand                 ! equation (A.25) in [1]

  end function Z_hc
  ! ############################################################################




  ! ############################################################################
  ! Function for the calculation of the hard-chain contribution to the
  !  chemical potential 'mu_hc(i)/kT' according to
  !  equation (A.35)  in [1].
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}
  !
  ! pure real function (in<PC-SAFT parameter>, in<densities>, in<molar fractions>, in<component index>, in<component index>)
  !  in  < saft_para :: saft_parameter {m(N_comp), d(N_comp)}>
  !  in  < rho       :: real(N_comp) >
  !  in  < x         :: real(N_comp) >
  ! ----------------------------------------------------------------------------
  pure real(dp) function beta_mu_hc(saft_para, rho, x)
    ! Input variables
    type(saft_parameter)                        :: saft_para          ! parameter container type
    real(dp), dimension(:), intent(in)          :: rho                ! density
    real(dp), intent(in)                        :: x                  ! molefraction

    ! Local variables
    real(dp)                                    :: m_dash             ! mean segment number
    real(dp), dimension(0:3)                    :: zeta               ! weighted densities
    real(dp), dimension(0:3, :), allocatable    :: zetax              ! partial densities
    real(dp), dimension(:),      allocatable    :: da_hs_dx           ! partial derivative of the reduced Helmholtz energy hard-sphere contribution w.r.t. the molar fractions
    real(dp)                                    :: summand = 0.0_dp   ! summand, help variable
    real(dp), dimension(:),      allocatable    :: da_hc_dx           ! partial derivative of the reduced Helmholtz energy hard-chain contribution w.r.t. the molar fractions

    integer                                     :: i, k               ! iteration variables


    ! Allocation
    allocate(da_hc_dx(saft_para%N_comp), zetax(0:3, saft_para%N_comp), da_hs_dx(saft_para%N_comp) )

    ! Mean segment number
    m_dash = sum(saft_para%m*x)                                     ! equation (A.5) in [1]

    ! Calculation of the summand according to equation (A.35) in [1]
    do i = 1, saft_para%N_comp
      summand = summand + x(i)*(m(i)-1._dp)/g_hs_ij(saft_para, rho, x, i, i) &
              & * dg_hs_ij_dxk(saft_para, rho, x, i, i)
    end do

    ! Weighted densities
    do n = 0, 3
      zeta(n) = PI/6._dp * rho * sum(x * saft_para%m * saft_para%d**n)    ! equation (9) or (A.8) in [1]
    end do

    ! Partial densities
    do n = 0, 3
      zetax(n,k) = PI/6._dp * rho * saft_para%m(k) * saft_para%d(k)**n    ! equation (A.34) in [1]
    end do

    ! Partial derivative of the residual reduced Helmholtz energy 'a_hs' w.r.t.
    !  the molar fraction 'x_k' according to equation (A.36) in [1]
    do k = 1, saft_para%N_comp
      da_hs_dx(k) = -zetax(0,k)/zeta(0) * a_tilde_hs(saft_para,rho,x) &               ! equation (A.36) in [1]
            & + 1._dp* ( (3._dp*( zetax(1,k)*zeta(2) + zeta(1)*zetax(2,k) )) &
            & /(1._dp - zeta(3)) &
            & + (3._dp*zeta(1)*zeta(2)*zetax(3,k))/(1._dp-zeta(3))**2 &
            & + (3._dp*zeta(2)*zetax(2,k))/(zeta(3)*(1._dp-zeta(3))**2) &
            & + (zeta(2)**3*zetax(3,k)*(3._dp*zeta(3)-1._dp)) &
            & / (zeta(3)**2 * (1._dp-zeta(3))**3 ) &
            & + ((3._dp*zeta(2)**2*zetax(2,k)*zeta(3) - 2._dp*zeta(2)**3*zetax(3,k)) &
            & / (zeta(3)**3) - zeta(0,k)) * log(1._dp-zeta(3))   &
            & + (zeta(0) - zeta(2)**3/zeta(3)**2) * zetax(3,k)/(1._dp-zeta(3)) )
    end do

    ! Calculation of the partial derivative of the hard-chain contribution to
    !  the reduced Helmholtz energy 'a_hc' w.r.t. the molar fractions
    do k = 1, saft_para%N_comp
      da_hc_dx(k) = m(k)*a_tilde_hs(saft_para, rho, x) + m_dash*da_hs_dx(k) &   ! equation (A.35) in [1]
                  & - summand
    end do

    ! Calculation of the hard chain contribution to the chemical potential
    do k = 1, saft_para%N_comp
      beta_mu_hc(k) = a_tilde_hc(saft_para, rho, x) + (Z_hc - 1._dp) + da_hc_dx(k)
          & - sum( x * da_hc_dx )
    end do

 end function beta_mu_hc
  ! ############################################################################




end module hard_chains_mod
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
