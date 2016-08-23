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


module parameter_mod
    use kinds_mod, only: dp
    implicit none
    private
    public :: PI

    real(dp), parameter :: PI = 4._dp * atan(1.0)
end module parameter_mod
