module thstrength_module
    use matinv_module
    implicit none
    contains
  
    subroutine THstrength(f, E1_layer, E2_layer, theta_layer, G12_layer, nju_layer, &
                          E1_tens_layer_strength, E2_tens_layer_strength, E1_compr_layer_strength, &
                          E2_compr_layer_strength, G12_layer_strength, strain_xy0, n)
      ! Inputs
      integer, intent(in) :: n
      real*8, intent(in) :: E1_layer(n), E2_layer(n), theta_layer(n), G12_layer(n)
      real*8, intent(in) :: nju_layer(n)
      real*8, intent(in) :: E1_tens_layer_strength(n), E2_tens_layer_strength(n)
      real*8, intent(in) :: E1_compr_layer_strength(n), E2_compr_layer_strength(n)
      real*8, intent(in) :: G12_layer_strength(n)
      real*8, intent(in) :: strain_xy0(3)
  
      ! Outputs
      real*8, intent(out) :: f(n)
  
      ! Local variables
      integer :: i
      real*8 :: c, s, E1, E2, G12, nju12, E1_strength, E2_strength
      real*8 :: S_ply(3, 3), Q(3, 3), T(3, 3)
      real*8 :: strain_xy(3), strain_12(3), sigma_12(3)
  
      ! Initialize f to zero
      f = 0.0d0
  
      ! Loop over each layer
      do i = 1, n
          strain_xy = strain_xy0
          c = cos(theta_layer(i))
          s = sin(theta_layer(i))
          E1 = E1_layer(i)
          E2 = E2_layer(i)
          G12 = G12_layer(i)
          nju12 = nju_layer(i)
  
          ! S_12 matrix (Gibson eq. 2.25)
          S_ply(1, 1) = 1.0d0 / E1
          S_ply(2, 2) = 1.0d0 / E2
          S_ply(1, 2) = -nju12 / E1
          S_ply(2, 1) = S_ply(1, 2)
          S_ply(3, 3) = 1.0d0 / G12
          S_ply(1, 3) = 0.0d0
          S_ply(2, 3) = 0.0d0
          S_ply(3, 1) = 0.0d0
          S_ply(3, 2) = 0.0d0
  
          ! Q_12 matrix
          Q = matinv3(S_ply)
  
          ! T matrix (transformation matrix)
          ! T matrix (transformation matrix)
          T = reshape([c**2, s**2, 2.0d0 * c * s, &
                       s**2, c**2, -2.0d0 * c * s, &
                      -c * s, c * s, c**2 - s**2], [3, 3])
  
          ! strain_12
          strain_xy(3) = strain_xy(3) / 2.0d0
          strain_12 = matmul(T, strain_xy)
          strain_12(3) = strain_12(3) * 2.0d0
  
          ! sigma_12 in the layer
          sigma_12 = matmul(Q, strain_12)
  
          if (sigma_12(1) > 0.0d0) then
              E1_strength = E1_tens_layer_strength(i)
          else
              E1_strength = E1_compr_layer_strength(i)
          end if
  
          if (sigma_12(2) > 0.0d0) then
              E2_strength = E2_tens_layer_strength(i)
          else
              E2_strength = E2_compr_layer_strength(i)
          end if
  
          ! Tsai-Hill equation (Gibson eq. 4.14)
          f(i) = sqrt((sigma_12(1) / E1_strength)**2 - sigma_12(1) * sigma_12(2) / E1_strength**2 + &
                      (sigma_12(2) / E2_strength)**2 + (sigma_12(3) / G12_layer_strength(i))**2)
      end do
  
    end subroutine THstrength
  
  end module thstrength_module
  