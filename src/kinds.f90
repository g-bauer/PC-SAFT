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
! MODULE kinds_mod
!  This module contains the double precision compiler independent kind
!  kind definition.
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module kinds_mod
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: dp

    integer, parameter :: dp = kind(real64)
end module kinds_mod
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
