module dispersion_mod
    use kinds_mod, only: dp
    use parameter_mod, only: PI
    use saft_parameter_mod, only: saft_parameter
    implicit none
    private
    public :: a_tilde_disp

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

    ! dispersive contribution to the helmholtz energy according to Equiation (A.10) of
    ! [1]
    pure real(dp) function a_tilde_disp(saft_para, rho, Temp, x)
        type(saft_parameter), intent(in) :: saft_para ! parameter container type
        real(dp), intent(in) :: rho ! density
        real(dp), intent(in) :: Temp ! temperature
        real(dp), dimension(:), intent(in) :: x ! molefraction

        integer :: i, j
        real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: sigma_ij
        real(dp), dimension(saft_para%N_comp,saft_para%N_comp) :: eps_k_ij
        real(dp) :: m_dash
        real(dp) :: m2_eps_sig3_dash
        real(dp) :: m2_eps2_sig3_dash
        real(dp) :: eta
        real(dp) :: I1, I2
        real(dp) :: C1

        associate(p => saft_para)

            ! combining rules for saft parameters
            sigma_ij = 0._dp
            eps_k_ij = 0._dp
            m2_eps_sig3_dash = 0._dp
            m2_eps2_sig3_dash = 0._dp
            do i = 1, p%N_comp
                do j = 1, p%N_comp
                    sigma_ij(j,i) = 0.5*(p%m(j) + p%m(i))
                    eps_k_ij(j,i) = sqrt(p%eps_k(j)*p%eps_k(i)) * (1._dp - p%k_ij(j,i))
                    m2_eps_sig3_dash = m2_eps_sig3_dash &
                        + x(j) * x(i) * p%m(j) * p%m(i) &
                        * eps_k_ij(j,i)/Temp * sigma_ij(j,i)
                    m2_eps2_sig3_dash = m2_eps_sig3_dash &
                        + x(j) * x(i) * p%m(j) * p%m(i) &
                        * (eps_k_ij(j,i)/Temp)**2 * sigma_ij(j,i)
                end do
            end do

            ! mean segment number in the mixture
            m_dash = sum(p%m*x)

            ! packing fraction eta
            eta = PI/6._dp * rho * sum(x * p%m * p%d**3)
        end associate

        ! integrals of the perturbation theory
        I1 = 0.0
        I2 = 0.0
        do i = 0, 6
            I1 = I1 + (a(0,i) &
            + (m_dash - 1._dp)/m_dash*a(1,i) &
            + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * a(2,i)) * eta**i
            I2 = I2 + (b(0,i) &
            + (m_dash - 1._dp)/m_dash*b(1,i) &
            + (m_dash - 1._dp)/m_dash * (m_dash - 2._dp)/m_dash * b(2,i)) * eta**i
        end do

        ! Compressibility expression C1
        C1 = 1._dp / (1._dp + m_dash*(8._dp*eta - 2._dp*eta**2)/(1._dp - eta)**4 &
        + (1._dp - m_dash)*(20._dp*eta - 27._dp*eta**2 + 12._dp*eta**3 &
        - 2._dp*eta**4)/((1._dp - eta)*(2._dp - eta))**2)

        a_tilde_disp = -2._dp*PI*rho*I1*m2_eps_sig3_dash &
                     - PI*rho*m_dash*C1*I2*m2_eps2_sig3_dash
    end function a_tilde_disp
end module dispersion_mod
