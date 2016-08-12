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
! MODULE saft_parameter_mod
!  This module contains the derived type for the PC-SAFT parameter and the
!  initialization subroutine:
!
!   - type saft_parameter
!     ~ <integer>                        :: N_comp       ! number of components in the system
!     ~ <real, dimension(N_comp)>        :: m            ! PC-SAFT parameter: segment number
!     ~ <real, dimension(N_comp)>        :: sigma        ! PC-SAFT parameter: segment diameter
!     ~ <real, dimension(N_comp)>        :: R            ! temperature dependent hard-sphere radius
!     ~ <real, dimension(N_comp)>        :: d            ! temperature dependent hard-sphere diameter
!     ~ <real, dimension(N_comp)>        :: eps_k        ! PC-SAFT parameter: energy parameter
!     ~ <real, dimension(N_comp,N_comp)> :: k_ij         ! adjustment for the binary interaction parameter
!
!
!   - subroutine init_saft_parameters(saft_para, Temp, m, sigma, eps_k, k_ij)
!     ~ <type(saft_parameter)>           :: saft_para    ! derived type
!     ~ <real>                           :: Temp         ! temperature
!     ~ <real, dimension(N_Comp)>        :: m            ! PC-SAFT parameter: segment number
!     ~ <real, dimension(N_comp)>        :: sigma        ! PC-SAFT parameter: segment diameter
!     ~ <real, dimension(N_comp)>        :: eps_k        ! PC-SAFT parameter: energy parameter
!     ~ <real, dimension(N_comp,N_comp)> :: k_ij         ! adjustment for the binary interaction parameter
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module saft_parameter_mod
  use kinds_mod, only: dp
  implicit none
  private
  public :: saft_parameter
  public :: init_saft_parameters
  public :: set_temperature_dependent_diameter


  ! parameter container type
  ! dimension = number of components
  type saft_parameter
    integer                               :: N_comp   ! number of components
    real(dp), allocatable, dimension(:)   :: m        ! PC-SAFT parameter: segment number
    real(dp), allocatable, dimension(:)   :: sigma    ! PC-SAFT parameter: segment diameter
    real(dp), allocatable, dimension(:)   :: R        ! segment radius (temperature dependent)
    real(dp), allocatable, dimension(:)   :: d        ! segment diameter (temperature dependent)
    real(dp), allocatable, dimension(:)   :: eps_k    ! PC-SAFT parameter: energy parameter eps/kB
    real(dp), allocatable, dimension(:,:) :: k_ij     ! adjustment for the binary interaction parameter
    real(dp) :: Temp
  end type saft_parameter

contains

  ! ############################################################################
  ! Subroutine for the initialization of the PC-SAFT parameter
  !
  ! subroutine (out<PC-SAFT parameter>, in<temperature>, in<segment numbers>, in<segment diameters>, in<energy parameter>, in<binary adjustment>)
  !  out < saft_para    :: saft_parameter      >
  !  in  < Temp         :: real                >
  !  in  < m            :: real(N_comp)        >
  !  in  < sigma        :: real(N_comp)        >
  !  in  < eps_k        :: real(N_comp)        >
  !  in  < k_ij         :: real(N_comp,N_Comp) >
  ! ----------------------------------------------------------------------------
  subroutine init_saft_parameters(saft_para, Temp, m, sigma, eps_k, k_ij)
    type(saft_parameter), intent(inout)   :: saft_para        ! derived type
    real(dp), intent(in)                  :: Temp             ! temperature
    real(dp), dimension(:), intent(in)    :: m                ! PC-SAFT parameter: segment numbers
    real(dp), dimension(:), intent(in)    :: sigma            ! PC-SAFT parameter: segment diameter
    real(dp), dimension(:), intent(in)    :: eps_k            ! PC-SAFT parameter: energy parameter
    real(dp), dimension(:,:), intent(in)  :: k_ij             ! adjustment to the binary interaction interaction parameter

    ! Get number of components from array size
    saft_para%N_comp = size(m)

    ! Allocate arrays
    allocate(saft_para%m(saft_para%N_comp))
    allocate(saft_para%sigma(saft_para%N_comp))
    allocate(saft_para%R(saft_para%N_comp))
    allocate(saft_para%eps_k(saft_para%N_comp))
    allocate(saft_para%k_ij(saft_para%N_comp,saft_para%N_comp))

    ! Copy values
    saft_para%m     = m
    saft_para%sigma = sigma
    saft_para%eps_k = eps_k
    saft_para%k_ij  = k_ij

    ! set temperature dependend values: d(T), R(T)
    call set_temperature_dependent_diameter(saft_para, Temp)
  end subroutine init_saft_parameters
  ! ############################################################################


  ! ############################################################################
  ! Subroutine for the initialization of temperature dependent segment
  !  diameter and radius.
  !
  ! subroutine (out<PC-SAFT parameter>, in<temperature>, in<segment numbers>, in<segment diameters>, in<energy parameter>, in<binary adjustment>)
  !  out < saft_para    :: saft_parameter{d(N_comp), R(N_comp)} >
  !  in  < Temp         :: real >
  ! ----------------------------------------------------------------------------
  subroutine set_temperature_dependent_diameter(saft_para, Temp)
    type(saft_parameter), intent(inout) :: saft_para
    real(dp), intent(in) :: Temp

    saft_para%d = saft_para%sigma* &
                & (1._dp - 0.12_dp*exp(-3._dp*saft_para%eps_k/Temp))
    saft_para%R = saft_para%d/2._dp
  end subroutine set_temperature_dependent_diameter
  ! ############################################################################


end module saft_parameter_mod
