module TWstrength_module
    use matinv_module
    implicit none
contains
    subroutine TWstrength(E1_layer, E2_layer, theta_layer, G12_layer, nju12_layer, F1t_layer, &
        F2t_layer, F1c_layer, F2c_layer, F12_layer, strain_xy0, FOS, n)
        
        implicit none

        ! Inputs
        integer, intent(in) :: n
        real*8, dimension(n), intent(in) :: E1_layer, E2_layer, theta_layer
        real*8, dimension(n), intent(in) :: G12_layer, nju12_layer, F1t_layer, F2t_layer
        real*8, dimension(n), intent(in) :: F1c_layer, F2c_layer, F12_layer
        real*8, dimension(3), intent(in) :: strain_xy0

        ! Output
        real*8, dimension(n), intent(out) :: FOS

        ! Local variables
        integer :: i
        real*8 :: c, s, E1, E2, G12, nju12
        real*8 :: S11, S12, S22, S66
        real*8, dimension(3, 3) :: S_matrix, Q, T
        real*8, dimension(3) :: strain_xy, strain_12, sigma_12
        real*8 :: f1, f2, f11, f22, f66, f12, a, b

        do i = 1, n
            strain_xy = strain_xy0
            c = cos(theta_layer(i))
            s = sin(theta_layer(i))
            E1 = E1_layer(i)
            E2 = E2_layer(i)
            G12 = G12_layer(i)
            nju12 = nju12_layer(i)

            ! S_12 matrix
            S11 = 1.0d0 / E1
            S12 = -nju12 / E1
            S22 = 1.0d0 / E2
            S66 = 1.0d0 / G12
            S_matrix = reshape([S11, S12, 0.0d0, S12, S22, 0.0d0, 0.0d0, 0.0d0, S66], shape=[3, 3])

            ! Q_12 matrix
            Q = matinv3(S_matrix)! Invert S to get Q (requires a matrix inversion subroutine)

            ! Transformation matrix T
            T = reshape([c**2, s**2, 2*c*s, &
                s**2, c**2, -2*c*s, &
                -c*s, c*s, c**2 - s**2], shape=[3, 3])

            ! Convert strain_xy to strain_12
            strain_xy(3) = strain_xy(3)/2.0d0
            strain_12 = matmul(T, strain_xy)
            strain_12(3) = 2.0 * strain_12(3) ! Adjust gamma for engineering strain

            ! Compute sigma_12
            sigma_12 = matmul(Q, strain_12)

            ! Compute f1, f2, f11, etc.
            f1 = 1.0 / F1t_layer(i) - 1.0 / F1c_layer(i)
            f2 = 1.0 / F2t_layer(i) - 1.0 / F2c_layer(i)
            f11 = 1.0 / (F1t_layer(i) * F1c_layer(i))
            f22 = 1.0 / (F2t_layer(i) * F2c_layer(i))
            f66 = 1.0 / (F12_layer(i)**2)
            if (sigma_12(1) > 0.0) then
                f12 = 1.0 / (2.0 * F1t_layer(i)) * (1.0 / F2c_layer(i) - 1.0 / F2t_layer(i))
            else
                f12 = -1.0 / (4.0 * F1c_layer(i) * F1t_layer(i))
            end if

            ! Compute a and b
            a = f11 * sigma_12(1)**2 + f22 * sigma_12(2)**2 + f66 * sigma_12(3)**2 + &
                2.0d0 * f12 * sigma_12(1) * sigma_12(2)
            b = f1 * sigma_12(1) + f2 * sigma_12(2)

            ! Compute Factor of Safety (FOS)
            FOS(i) = (-b + sqrt(b**2 + 4.0d0 * a)) / (2.0d0 * a)
        end do
    end subroutine
end module TWstrength_module
