module laminate_module
  use matinv_module
  implicit none
contains

  subroutine laminate(A_lam, D_lam, E1_layer, E2_layer, theta_layer, t_layer, G12_layer, nju_layer, z0, n)
    ! Inputs
    integer, intent(in) :: n
    real*8, intent(in) :: E1_layer(n), E2_layer(n), theta_layer(n), t_layer(n)
    real*8, intent(in) :: G12_layer(n), nju_layer(n)
    real*8, intent(in) :: z0

    ! Outputs
    real*8, intent(out) :: A_lam(3, 3), D_lam(3, 3)

    ! Local variables
    integer :: i
    real*8 :: c, s, E1, E2, G12, nju12
    real*8 :: S_ply(3, 3), Q(3, 3), T(3, 3), Q0(3, 3)
    real*8 :: S11, S22, S12, S21, S66
    real*8 :: zk, zk_1

    ! Initialize A_lam and D_lam matrices to zero
    A_lam = 0.0d0
    D_lam = 0.0d0

    ! Loop over each layer
    do i = 1, n
      c = cos(theta_layer(i))
      s = sin(theta_layer(i))
      E1 = E1_layer(i)
      E2 = E2_layer(i)
      G12 = G12_layer(i)
      nju12 = nju_layer(i)

      ! S_12 matrix (Gibson eq. 2.25)
      S11 = 1.0d0 / E1
      S22 = 1.0d0 / E2
      S12 = -nju12 / E1
      S21 = S12
      S66 = 1.0d0 / G12
      S_ply = reshape([S11, S12, 0.0d0, S21, S22, 0.0d0, 0.0d0, 0.0d0, S66], [3, 3])

      ! Q_12 matrix
      Q = matinv3(S_ply)
      Q(:, 3) = Q(:, 3) * 2.0d0

      ! T matrix (transformation matrix)
      T = reshape([c**2, s**2, 2.0d0 * c * s, s**2, c**2, -2.0d0 * c * s, -c * s, c * s, c**2 - s**2], [3, 3])

      ! Q_xy matrix
      Q0 = matmul(matmul(matinv3(T), Q), T)
      Q0(:, 3) = Q0(:, 3) / 2.0d0

      ! A_xy of composite laminate (Gibson eq. 7.38)
      A_lam = A_lam + Q0 * t_layer(i)

      ! D_xy of composite laminate (Gibson eq. 7.40)
      zk = -sum(t_layer) + sum(t_layer(1:n-i+1)) + z0
      zk_1 = -sum(t_layer) + sum(t_layer(1:n-i)) + z0
      D_lam = D_lam + Q0 * (zk**3 - zk_1**3) / 3.0d0
    end do

  end subroutine laminate
  
end module laminate_module
  