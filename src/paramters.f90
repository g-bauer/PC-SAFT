module parameter_mod
    use kinds_mod, only: dp
    implicit none
    private
    public :: PI

    real(dp), parameter :: PI = 3.14159265359
end module parameter_mod
