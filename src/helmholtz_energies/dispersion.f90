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
! MODULE dispersion_mod
!  This module contains all dispersion contributions for the PC-SAFT equation
!  of state:
!
!   - pure real function a_tilde_disp(saft_para, rho, Temp, x)
!     ~ <type(saft_parameter)>    :: saft_para      ! PC-SAFT parameters
!     ~ <real, dimension(N_comp)> :: rho            ! densities
!     ~ <real>                    :: Temp           ! temperature
!     ~ <real, dimension(N_comp)> :: x              ! mole fractions
!
!   - subroutine ....
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module dispersion_mod
  use kinds_mod, only: dp                             ! loading 'real64' precision specifier (compiler independent)
  use parameter_mod, only: PI                         ! loading pi=3.14159
  use saft_parameter_mod, only: saft_parameter        ! loading derived type 'saft_parameter'
  implicit none
  private                                             ! declaring everythin in this module 'private'
  public :: a_tilde_disp                              ! specifically declaring 'public'

  ! Universal model constants for integrals of perturbation for hard chains.
  !  Table 1, used for equations (18) & (19) in
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}
  real(dp), parameter, dimension(0:2,0:6) :: a = &
    reshape([   0.9105631445, -0.3084016918, -0.0906148351, &
    0.6361281449,  0.1860531159,  0.4527842806, &
    2.686134789,  -2.503004726,   0.5962700728, &
    -26.54736249,   21.41979363,   -1.724182913, &
    97.75920878,  -65.25588533,   -4.130211253, &
    -159.5915409,    83.31868048,   13.77663187, &
    91.29777408,  -33.74692293,   -8.672847037 ],[3,7])
    real(dp), parameter, dimension(0:2,0:6) :: b = &
    reshape([   0.7240946941, -0.5755498075,  0.0976883116, &
    2.238279186,   0.6995095521, -0.2557574982, &
    -4.002584949,   3.892567339,  -9.155856153, &
    -21.00357682,  -17.21547165,   20.64207597, &
    26.85564136,  192.6722645,   -38.80443005, &
    206.5513384,  -161.8264617,    93.62677408, &
    -355.6023561,  -165.2076935,   -29.66690559 ],[3,7])

contains
  ! ############################################################################
  ! Function for the calculation of the dispersive contribution to the reduced
  !  Helmholtz energy 'a_tilde_disp = A_disp/NkT' according to equation (A.10) in
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}

  ! pure real function (in<PC-SAFT parameter>, in<densites>, in<temperature>, in<mole fractions>)
  !  in  < saft_para    :: saft_parameter {m(N_comp), sigma(N_comp), R(N_Comp),
  !                         d(N_comp), eps_k(N_comp), k_ij(N_Comp,N_comp))} >
  !  in  < rho          :: real(N_comp) >
  !  in  < Temp         :: real         >
  !  in  < x            :: real(N_comp) >
  !  out < a_tilde_disp :: real         >
  ! ----------------------------------------------------------------------------
  pure real(dp) function a_tilde_disp(saft_para, rho, Temp, x)
    ! Input variables
    type(saft_parameter),   intent(in) :: saft_para     ! parameter container type
    real(dp),               intent(in) :: rho           ! density
    real(dp),               intent(in) :: Temp          ! temperature
    real(dp), dimension(:), intent(in) :: x             ! molefraction

    ! Local variables
    integer :: i, j                                     ! iteration variables

    real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: sigma_ij ! array for mean PC-SAFT diameter
    real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: eps_k_ij ! array for binary energy interactions
    real(dp) :: m_dash                                  ! mean segment number
    real(dp) :: m2_eps_sig3_dash                        ! mean (m(i)^2 eps_k(i)/T sig(i)^3)
    real(dp) :: m2_eps2_sig3_dash                       ! mean (m(i)^2 (eps_k(i)/T)^2 sig(i)^3)
    real(dp) :: eta                                     ! packing fraction
    real(dp) :: I1, I2                                  ! integrals for chain perturbations; equations (14), (15) in [1]
    real(dp) :: C1                                      ! derivative of compressibilities; equations (13) in [1]


    associate(p => saft_para)

      ! Combination rules for PC-SAFT parameters
      sigma_ij = 0._dp                ! zero initialization for summation
      eps_k_ij = 0._dp                ! zero initialization for summation
      m2_eps_sig3_dash = 0._dp        ! zero initialization for summation
      m2_eps2_sig3_dash = 0._dp       ! zero initialization for summation

      ! Calculation of combination rules
      do i = 1, p%N_comp
        do j = 1, p%N_comp
          sigma_ij(j,i) = 0.5*(p%m(j) + p%m(i))                       ! mean segment diameter; equation (23) or (A.14) in [1]
          eps_k_ij(j,i) = sqrt(p%eps_k(j)*p%eps_k(i)) &
                        & * (1._dp - p%k_ij(j,i))                     ! binary interaction parameter; equation (24) or (A.15) in [1]
          m2_eps_sig3_dash = m2_eps_sig3_dash &                       ! abbreviation (A.12) in [1]
                           & + x(j) * x(i) * p%m(j) * p%m(i) &
                           & * eps_k_ij(j,i)/Temp * sigma_ij(j,i)**3
          m2_eps2_sig3_dash = m2_eps_sig3_dash &                      ! abbreviation (A.13) in [1]
                            & + x(j) * x(i) * p%m(j) * p%m(i) &
                            & * (eps_k_ij(j,i)/Temp)**2 * sigma_ij(j,i)**3
        end do
      end do

      ! Calculation of mean segment number in the mixture
      m_dash = sum(p%m*x)                                 ! equation (6) or (A.5) in [1]

      ! Calculation of packing fraction eta
      eta = PI/6._dp * rho * sum(x * p%m * p%d**3)        ! from equation (A.20) in [1]
    end associate

    ! Calculation of the integrals of the perturbation theory for hard chains
    I1 = 0.0_dp        ! zero initialization for summation
    I2 = 0.0_dp        ! zero initialization for summation
    do i = 0, 6
      I1 = I1 + (a(0,i) &                                                     ! combinations of equations (16) & (18) or equations (A.16) & (A.18) in [1]
         & + (m_dash - 1._dp)/m_dash*a(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * a(2,i)) * eta**i
      I2 = I2 + (b(0,i) &                                                     ! combinations of equations (17) & (19) or equations (A.17) & (A.19) in [1]
         & + (m_dash - 1._dp)/m_dash*b(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * b(2,i)) * eta**i
    end do

    ! Calculation of the compressibility expression C1
    C1 = 1._dp / (1._dp + m_dash*(8._dp*eta - 2._dp*eta**2)/(1._dp - eta)**4 &  ! equation (13), (25) or (A.11); attention: equation (A.11) is missing the reziprocal (-)^-1 on the last line
       & + (1._dp - m_dash)*(20._dp*eta - 27._dp*eta**2 + 12._dp*eta**3 &
       & - 2._dp*eta**4)/((1._dp - eta)*(2._dp - eta))**2)

    ! Calculation of the reduced Helmholtz energy
    a_tilde_disp = -2._dp*PI*rho*I1*m2_eps_sig3_dash &            ! equation (10) oder (A.10) in [1]
                 & - PI*rho*m_dash*C1*I2*m2_eps2_sig3_dash
  end function a_tilde_disp
    ! ##########################################################################



  ! ############################################################################
  ! Function for the calculation of the dispersion compressibility 'Z_hs'
  !  according to equations (A.28) in [1].
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}
  !
  ! pure real function (in<PC-SAFT parameter>, in<densities>, in<temperature>, in<molar fractions>)
  !  in  < saft_para :: saft_parameter {m(N_comp), d(N_comp)}>
  !  in  < rho       :: real(N_comp) >
  !  in  < Temp      :: real         >
  !  in  < x         :: real(N_comp) >
  ! ----------------------------------------------------------------------------
  pure real(dp) function Z_disp(saft_para, rho, Temp, x)
    ! Input variables
    type(saft_parameter), intent(in)   :: saft_para     ! parameter container type
    real(dp), intent(in)               :: rho           ! density
    real(dp), intent(in)               :: Temp          ! temperature
    real(dp), dimension(:), intent(in) :: x             ! molefraction

    ! Local variables
    integer :: i, j                                     ! iteration variables

    real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: sigma_ij ! array for mean PC-SAFT diameter
    real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: eps_k_ij ! array for binary energy interactions
    real(dp) :: m_dash                                  ! mean segment number
    real(dp) :: m2_eps_sig3_dash                        ! mean (m(i)^2 eps_k(i)/T sig(i)^3)
    real(dp) :: m2_eps2_sig3_dash                       ! mean (m(i)^2 (eps_k(i)/T)^2 sig(i)^3)
    real(dp) :: eta                                     ! packing fraction
    real(dp) :: I1, I2                                  ! integrals for chain perturbations; equations (14), (15) in [1]
    real(dp) :: C1, C2                                  ! derivative of compressibilities; equations (13) in [1]


    associate(p => saft_para)
      ! Combination rules for PC-SAFT parameters
      sigma_ij = 0._dp                ! zero initialization for summation
      eps_k_ij = 0._dp                ! zero initialization for summation
      m2_eps_sig3_dash = 0._dp        ! zero initialization for summation
      m2_eps2_sig3_dash = 0._dp       ! zero initialization for summation

      ! Calculation of combination rules
      do i = 1, p%N_comp
        do j = 1, p%N_comp
          sigma_ij(j,i) = 0.5*(p%m(j) + p%m(i))                       ! mean segment diameter; equation (23) or (A.14) in [1]
          eps_k_ij(j,i) = sqrt(p%eps_k(j)*p%eps_k(i)) &
                        & * (1._dp - p%k_ij(j,i))                     ! binary interaction parameter; equation (24) or (A.15) in [1]
          m2_eps_sig3_dash = m2_eps_sig3_dash &                       ! abbreviation (A.12) in [1]
                           & + x(j) * x(i) * p%m(j) * p%m(i) &
                           & * eps_k_ij(j,i)/Temp * sigma_ij(j,i)**3
          m2_eps2_sig3_dash = m2_eps_sig3_dash &                      ! abbreviation (A.13) in [1]
                            & + x(j) * x(i) * p%m(j) * p%m(i) &
                            & * (eps_k_ij(j,i)/Temp)**2 * sigma_ij(j,i)**3
        end do
      end do

      ! Calculation of mean segment number in the mixture
      m_dash = sum(p%m*x)                                 ! equation (6) or (A.5) in [1]

      ! Calculation of packing fraction eta
      eta = PI/6._dp * rho * sum(x * p%m * p%d**3)        ! from equation (A.20) in [1]
    end associate

    ! Calculation of the integrals of the perturbation theory for hard chains
    !  and their derivatives w.r.t eta
    I2     = 0.0        ! zero initialization for summation
    detaI1 = 0.0        ! zero initialization
    detaI2 = 0.0        ! zero initialization
    do i = 0, 6
      I2 = I2 + (b(0,i) &                                                     ! combinations of equations (17) & (19) or equations (A.17) & (A.19) in [1]
         & + (m_dash - 1._dp)/m_dash*b(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * b(2,i)) * eta**i
      detaI1 = detaI1 + (a(0,i)                 &                             ! equation (A.29) in [1]
             & + (m_dash - 1._dp)/m_dash*a(1,i) &
             & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * a(2,i)) &
             & * (real(i) + 1._dp) * eta**i
      detaI2 = detaI2 +(b(0,i)                 &                              ! equation (A.30) in [1]
             & + (m_dash - 1._dp)/m_dash*b(1,i) &
             & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * b(2,i)) &
             & * (real(i) + 1._dp) * eta**i
    end do

    C1 = 1._dp / (1._dp + m_dash*(8._dp*eta - 2._dp*eta**2)/(1._dp - eta)**4 &  ! equation (13), (25) or (A.11); attention: equation (A.11) is missing the reziprocal (-)^-1 on the last line
       & + (1._dp - m_dash)*(20._dp*eta - 27._dp*eta**2 + 12._dp*eta**3 &
       & - 2._dp*eta**4)/((1._dp - eta)*(2._dp - eta))**2)
    C2 = - C1**2 * ( m_dash*(-4._dp*eta**2 + 20._dp*eta + 8._dp)/(1._dp-eta)**5  & ! equation (A.31) in [1]
       & (1._dp - m_dash) * (2._dp*eta**3 + 12._dp*eta**2 - 48._dp*eta + 40._dp) &
       &  / ((1._dp-eta)*(2._dp-eta))**3 )

    ! Calculation of the dispersion compressibility contribution
    Z_disp = -2._dp*PI*rho*detaI1*m2_eps_sig3_dash &
           & -PI*rho*m_dash* ( C1*detaI2 + C2*eta*I2 ) * m2_eps2_sig3_dash

  end function Z_disp


  ! ############################################################################
  ! Function for the calculation of the dispersive contribution to the
  !  chemical potential 'mu_hc(i)/kT' according to
  !  equation (A.33) & (A.38)  in [1].
  !  {Gross, J., Sadowski, G.: Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules
  !  Industrial Engineering & Chemistry Research, 2001, Vol. 40 (4), 1244-1260}
  !
  ! pure real function (in<PC-SAFT parameter>, in<densities>, in<molar fractions>, in<component index>, in<component index>)
  !  in  < saft_para :: saft_parameter {m(N_comp), d(N_comp)}>
  !  in  < rho       :: real(N_comp) >
  !  in  < x         :: real(N_comp) >
  ! ----------------------------------------------------------------------------
  pure real(dp) function beta_mu_hc(saft_para, rho, Temp, x)
    ! Input variables
    type(saft_parameter)                        :: saft_para      ! parameter container type
    real(dp), dimension(:), intent(in)          :: rho            ! density
    real(dp), intent(in)                        :: Temp           ! temperature
    real(dp), intent(in)                        :: x              ! molefraction

    ! Local variables
    real(dp)                                    :: eta            ! packing fractions
    real(dp), dimension(:), allocatable         :: zeta3          ! partial derivative of packing fraction w.r.t. molar fraction
    real(dp)                                    :: I1, I2         ! integrals for chain perturbations; equations (14), (15) in [1]
    real(dp)                                    :: C1, C2         ! derivative of compressibilities; equations (13) in [1]
    real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: sigma_ij ! array for mean PC-SAFT diameter
    real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: eps_k_ij ! array for binary energy interactions
    real(dp)                                    :: m_dash         ! mean segment number
    real(dp), dimension(:), allocatable         :: I1x            ! derivative of 'I1' w.r.t. molar fraction
    real(dp), dimension(:), allocatable         :: I2x            ! derivative of 'I2' w.r.t. molar fraction
    real(dp), dimension(:), allocatable         :: C1x            ! derivative of 'C1' w.r.t. molar fraction
    real(dp)                                    :: m2_eps2_sig3_dash      ! mean (m(i)^2 eps_k(i)/T sig(i)^3)
    real(dp)                                    :: m2_eps2_sig3_dash      ! mean (m(i)^2 (eps_k(i)/T)^2 sig(i)^3)
    real(dp), dimension(:), allocatable         :: m2_eps_sig3_dash_x     ! derivative of 'm2_eps2_sig3_dash' w.r.t. molar fraction
    real(dp), dimension(:), allocatable         :: m2_eps2_sig3_dash_x    ! derivative of 'm2_eps2_sig3_dash' w.r.t. molar fraction
    real(dp), dimension(:,:), allocatable       :: ax             ! derivative of 'a' w.r.t. molar fraction
    real(dp), dimension(:,:), allocatable       :: bx             ! derivative of 'b' w.r.t. molar fraction
    real(dp), dimension(:), allocatable         :: da_disp_dx     ! derivative of 'a_tilde_disp' w.r.t. molar fractions

    integer                                     :: i, k           ! iteration variables


    ! Allocation
    allocate( zeta3(N_comp), I1x(saft_para%N_comp), I2x(saft_para%N_comp), &
            & C1x(saft_para%N_comp), &
            & m2_eps_sig3_dash_x(saft_para%N_comp), &
            & m2_eps2_sig3_dash_x(saft_para%N_comp), &
            & ax(0:6, saft_para%N_comp), bx(0:6, saft_para%N_comp) )

    ! Mean segment number
    m_dash = sum(saft_para%m*x)                                     ! equation (A.5) in [1]

    ! Calculation of combination rules
    do i = 1, p%N_comp
      do j = 1, p%N_comp
        sigma_ij(j,i) = 0.5*(p%m(j) + p%m(i))                       ! mean segment diameter; equation (23) or (A.14) in [1]
        eps_k_ij(j,i) = sqrt(p%eps_k(j)*p%eps_k(i)) &
                      & * (1._dp - p%k_ij(j,i))                     ! binary interaction parameter; equation (24) or (A.15) in [1]
        m2_eps_sig3_dash = m2_eps_sig3_dash &                       ! abbreviation (A.12) in [1]
                         & + x(j) * x(i) * p%m(j) * p%m(i) &
                         & * eps_k_ij(j,i)/Temp * sigma_ij(j,i)**3
        m2_eps2_sig3_dash = m2_eps_sig3_dash &                      ! abbreviation (A.13) in [1]
                          & + x(j) * x(i) * p%m(j) * p%m(i) &
                          & * (eps_k_ij(j,i)/Temp)**2 * sigma_ij(j,i)**3
      end do
    end do

    ! Calculation of packing fraction 'eta' and 'zeta3'
    eta   = PI/6._dp * rho * sum(x * p%m * p%d**3)                  ! equation (A.20) in [1]
    zeta3 = PI/6._dp * rho * saft_para%m * saft_para%d**3           ! equation (A.34) in [1]

    ! Calculation of the integrals of the perturbation theory for hard chains
    !  and their derivatives w.r.t eta
    I1     = 0.0_dp     ! zero initialization for summation
    I2     = 0.0_dp     ! zero initialization for summation
    detaI1 = 0.0_dp     ! zero initialization for summation
    detaI2 = 0.0_dp     ! zero initialization for summation
    do i = 0, 6
      I1 = I1 + (a(0,i) &                                                       ! combinations of equations (16) & (18) or equations (A.16) & (A.18) in [1]
         & + (m_dash - 1._dp)/m_dash*a(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * a(2,i)) * eta**i
      I2 = I2 + (b(0,i) &                                                       ! combinations of equations (17) & (19) or equations (A.17) & (A.19) in [1]
         & + (m_dash - 1._dp)/m_dash*b(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * b(2,i)) * eta**i
      ax(i,:) = m / m_dash * a(1,i) + m / m_dash * a(2,i)                       ! equation (A.44) in [1]
      bx(i,:) = m / m_dash * b(1,i) + m / m_dash * b(2,i)                       ! equation (A.45) in [1]
      I1x = I1x + (a(0,i) &                                                     ! equation (A.42) in [1]
         & + (m_dash - 1._dp)/m_dash*a(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * a(2,i)) &
         & * i * zeta3 * eta**(i-1) &
         & + ax(i,:) * eta**i
      I2x = I2x + (b(0,i) &                                                     ! equation (A.43) in [1]
         & + (m_dash - 1._dp)/m_dash*b(1,i) &
         & + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * b(2,i)) &
         & * i * zeta3 * eta**(i-1) &
         & + ax(i,:) * eta**i
    end do

    ! Calculation of C's
    C1 = 1._dp / (1._dp + m_dash*(8._dp*eta - 2._dp*eta**2)/(1._dp - eta)**4 &  ! equation (13), (25) or (A.11); attention: equation (A.11) is missing the reziprocal (-)^-1 on the last line
       & + (1._dp - m_dash)*(20._dp*eta - 27._dp*eta**2 + 12._dp*eta**3 &
       & - 2._dp*eta**4)/((1._dp - eta)*(2._dp - eta))**2)
    C2 = - C1**2 * ( m_dash*(-4._dp*eta**2 + 20._dp*eta + 8._dp)/(1._dp-eta)**5  & ! equation (A.31) in [1]
       & (1._dp - m_dash) * (2._dp*eta**3 + 12._dp*eta**2 - 48._dp*eta + 40._dp) &
       &  / ((1._dp-eta)*(2._dp-eta))**3 )
    C1x = C2*zeta3 - C1**2* ( m*(8._dp*eta - 2._dp*eta**2)/(1._dp-eta)**4  &
        & - m*(20._dp*eta-27._dp*eta**2+12._dp*eta**3-2._dp*eta**4) &
        & / ((1._dp-eta)*(2._dp-eta))**2 )

    ! Calculation of abbreviations
    do k = 1, saft_para%N_Comp
      m2_eps_sig3_dash_x(k)  = 2._dp*saft_para%m(k) &
                             & * sum( x*m*(eps_k_ij(:,k)/Temp)   *sigma_ij(:,k) )
      m2_eps2_sig3_dash_x(k) = 2._dp*saft_para%m(k) &
                             & * sum( x*m*(eps_k_ij(:,k)/Temp)**2*sigma_ij(:,k) )
    end do

    ! Calculation of the derivative of the dispersion contribution to the
    !  reduced Helmholtz energy w.r.t. the molar fractions
    da_disp_dx = -2._dp*PI* (I1x*m2_eps_sig3_dash + I1*m2_eps_sig3_dash_x ) &           ! equation (A.38) in [1]
               & - PI*rho*( (m*C1*I2+m_dash*C1x*I2+m_dash*C1*I2x) * m2_eps2_sig3_dash &
               & + m_dash*C1*I2*m2_eps2_sig3_dash_x )

    ! Chemical potential contribution of the dispersive interaction
    do k = 1, saft_para%N_comp
      beta_mu_hs(k) = a_tilde_disp(saft_para, rho, Temp, x) + (Z_disp - 1._dp) &          ! equation (A.33) in [1]
                    & + da_disp_dx(k) - sum( x * da_disp_dx )
    end do


 end function beta_mu_hc
  ! ############################################################################

end module dispersion_mod
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
