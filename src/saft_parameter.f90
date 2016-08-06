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
        integer :: N_comp
        real(dp), allocatable, dimension(:) :: m ! segment number
        real(dp), allocatable, dimension(:) :: sigma ! segment diameter
        real(dp), allocatable, dimension(:) :: R ! segment radius (temperature dependent)
        real(dp), allocatable, dimension(:) :: d ! segment diameter (temperature dependent)
        real(dp), allocatable, dimension(:) :: eps_k ! energy parameter eps/kB
        real(dp), allocatable, dimension(:,:) :: k_ij ! binary parameter
        real(dp) :: Temp
    end type saft_parameter

contains

    ! allocate parameter arrays and set initial values
    subroutine init_saft_parameters(saft_para, Temp, m, sigma, eps_k, k_ij)
        type(saft_parameter), intent(inout) :: saft_para
        real(dp), intent(in) :: Temp
        real(dp), dimension(:), intent(in)  :: m
        real(dp), dimension(:), intent(in)  :: sigma
        real(dp), dimension(:), intent(in)  :: eps_k
        real(dp), dimension(:,:), intent(in)  :: k_ij

        ! get number of components from array size
        saft_para%N_comp = size(m)

        ! allocate arrays
        allocate(saft_para%m(saft_para%N_comp))
        allocate(saft_para%sigma(saft_para%N_comp))
        allocate(saft_para%R(saft_para%N_comp))
        allocate(saft_para%eps_k(saft_para%N_comp))
        allocate(saft_para%k_ij(saft_para%N_comp,saft_para%N_comp))

        ! copy values
        saft_para%m = m
        saft_para%sigma = sigma
        saft_para%eps_k = eps_k
        saft_para%k_ij = k_ij

        ! set temperature dependend values: d(T), R(T)
        call set_temperature_dependent_diameter(saft_para, Temp)
    end subroutine init_saft_parameters

    subroutine set_temperature_dependent_diameter(saft_para, Temp)
        type(saft_parameter), intent(inout) :: saft_para
        real(dp), intent(in) :: Temp

        saft_para%d = saft_para%sigma* &
        & (1._dp - 0.12_dp*exp(-3._dp*saft_para%eps_k/Temp))
        saft_para%R = saft_para%d/2._dp
    end subroutine set_temperature_dependent_diameter

end module saft_parameter_mod
