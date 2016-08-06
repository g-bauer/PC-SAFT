module kinds_mod
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: dp

    integer, parameter :: dp = kind(real64)
end module kinds_mod
