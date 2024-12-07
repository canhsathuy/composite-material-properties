program beam_without_core
    use ply_module
    use laminate_module
    use TWstrength_module
    use matinv_module
    implicit none
    
    real*8 :: E1c_282, E2c_282, G12c_282
    real*8 :: nju12c_282
    real*8 :: rho_282
    real*8 :: sLc_tens_282, sLc_compr_282
    real*8 :: sTc_tens_282, sTc_compr_282, sLTc_282
    real*8 :: m_282, Vf_282, V_282, tc_282, t_woven

    real*8 :: E1_woven, E2_woven, G12_woven
    real*8 :: nju12_woven

    ! AS4 Carbon fiber
    real*8 :: E1f_carb = 235d9
    real*8 :: E2f_carb = 15d9
    real*8 :: G12f_carb = 27d9
    real*8 :: nju12f_carb = 0.2d0
    real*8 :: sLf_tens_carb = 3.7d9
    !real*8 :: sLf_compr_carb = 3.3d9  ! (from T300 carbon fiber)
    real*8 :: rho_f_carb = 1810
    real*8 :: vf_carb = 0.63d0
    integer :: style_282 = 1  ! uniform distribution (1); layered distribution (0)

    ! Matrix (LY 564 / Aradur 560)
    real*8 :: Em = 4.3d9
    real*8 :: nju12_m = 0.35d0
    real*8 :: G12m
    real*8 :: sLm_tens = 69d6
    real*8 :: sLm_compr = 200d6  ! From Epoxy 3501-6
    real*8 :: sLTm = 100d6  ! From Epoxy 3501-6
    real*8 :: rho_m = 1270

    ! Woven lamina moduli
    real*8, dimension(2) :: E1_layer_woven, E2_layer_woven, G12_layer_woven
    real*8, dimension(2) :: theta_layer_woven, t_layer_woven, nju12_layer_woven
    real*8, dimension(2) :: F1t_layer_woven, F2t_layer_woven
    real*8, dimension(2) :: F1c_layer_woven, F2c_layer_woven, F12_layer_woven

    ! Unidirectional composite
    real*8 :: E1_ud, E2_ud, G12_ud
    real*8 :: nju12_ud
    real*8 :: rho_ud
    real*8 :: sL_tens_ud, sL_compr_ud
    real*8 :: sT_tens_ud, sT_compr_ud, sLT_ud
    real*8 :: m_ud, Vf_ud, V_ud, t_ud
    integer :: style_ud

    ! Khai bao bien cua spar
    real*8 :: h_spar, L, p
    integer :: n_layer_spar, n_layer_cap
    real*8, dimension(10) :: E1_layer_spar, E2_layer_spar, G12_layer_spar
    real*8, dimension(10) :: theta_layer_spar, nju12_layer_spar, rho_spar
    real*8, dimension(10) :: t_layer_spar
    real*8 :: F1t_spar, F1c_spar, F2t_spar, F2c_spar, F12_spar
    ! Khai bao bien cua cap
    real*8, dimension(16) :: E1_layer_cap, E2_layer_cap, G12_layer_cap
    real*8, dimension(16) :: theta_layer_cap, nju12_layer_cap, rho_cap
    real*8, dimension(16) :: t_layer_cap

    ! Intermediate and output variables
    real*8 :: t_cap, t_spar, z0, Ex_spar
    real*8, dimension(3) :: stress, strain
    real*8, dimension(2) :: FOS
    real*8 :: I_spar, EI_cap, EI, w_cap, w, d_max
    real*8, dimension(3,3) :: A_lam_woven, D_lam_woven, Q_lam_woven, S_lam_woven
    real*8, dimension(3,3) :: A_lam_cap, D_lam_cap, Q_lam_cap
    real*8, dimension(3,3) :: A_lam_spar, D_lam_spar, Q_lam_spar, S_lam_spar
    real*8, dimension(3) :: kapa_unit

    real*8, parameter :: pi = 3.14159265358979d0
    integer :: i

    real*8, dimension(101) :: x, M, kapa, sigma_cap

    ! Set material and geometry parameters
    n_layer_spar = 10
    n_layer_cap = 16
    h_spar = 169.053e-3
    ! core_thickness = 50e-3
    L = 8.0     ! Beam length
    p = 0.01d6  ! Pressure in Pa

    G12m = Em/2.0d0/(1 + nju12_m)
    call ply(E1c_282, E2c_282, G12c_282, nju12c_282, sLc_tens_282, sLc_compr_282, &
             sTc_tens_282, sTc_compr_282, sLTc_282, rho_282, E1f_carb, E2f_carb, &
             G12f_carb, Em, G12m, 0.62d0, sLf_tens_carb, sLm_tens, &
             sLm_compr, sLTm, rho_f_carb, rho_m, nju12f_carb, nju12_m, style_282)
    
    m_282 = 0.197 / 2.0         ! mass of fiber per m^2 in each ply
    Vf_282 = m_282 / rho_f_carb ! volume of fiber per m^2
    V_282 = Vf_282 / vf_carb    ! total volume of the ply per m^2
    tc_282 = V_282              ! ply thickness
    ! Module cua tung lop vai dan
    E1_layer_woven = (/ E1c_282, E1c_282 /)
    E2_layer_woven = (/ E2c_282, E2c_282 /)
    G12_layer_woven = (/ G12c_282, G12c_282 /)
    theta_layer_woven = (/ 0.0d0, pi/2.0 /)
    t_layer_woven = (/ tc_282, tc_282 /)
    nju12_layer_woven = (/ nju12c_282, nju12c_282 /)
    ! Do ben cua tung lop vai dan
    F1t_layer_woven = (/sLc_tens_282, sLc_tens_282/)
    F2t_layer_woven = (/sTc_tens_282, sTc_tens_282/)
    F1c_layer_woven = (/sLc_compr_282, sLc_compr_282/)
    F2c_layer_woven = (/sTc_compr_282, sTc_compr_282/)
    F12_layer_woven = (/sLTc_282, sLTc_282/)
    ! Tinh thong so vai dan
    call laminate(A_lam_woven, D_lam_woven, E1_layer_woven, E2_layer_woven, theta_layer_woven, &
        t_layer_woven, G12_layer_woven, nju12_layer_woven, 0.0d0, 2)
    
    Q_lam_woven = A_lam_woven / sum(t_layer_woven)
    S_lam_woven = matinv3(Q_lam_woven)
    
    E1_woven = 1.0 / S_lam_woven(1,1)  ! Gibson Eq. 2.25
    E2_woven = 1.0 / S_lam_woven(2,2)
    G12_woven = 1.0 / S_lam_woven(3,3)
    nju12_woven = -S_lam_woven(2,1) / S_lam_woven(1,1)
    t_woven = 2.0 * tc_282

    ! Tinh do ben cua vai dan
    ! Longitudinal Tensile strength
    stress = [1.0, 0.0, 0.0]
    strain = matmul(Q_lam_woven, stress)
    call TWstrength(E1_layer_woven, E2_layer_woven, theta_layer_woven, G12_layer_woven, nju12_layer_woven, &
                F1t_layer_woven, F2t_layer_woven, F1c_layer_woven, F2c_layer_woven, F12_layer_woven, &
                strain, FOS, 2)
    F1t_spar = maxval(FOS)

    ! Longitudinal Compressive strength
    stress = [-1.0, 0.0, 0.0]
    strain = matmul(Q_lam_woven, stress)
    call TWstrength(E1_layer_woven, E2_layer_woven, theta_layer_woven, G12_layer_woven, nju12_layer_woven, &
                F1t_layer_woven, F2t_layer_woven, F1c_layer_woven, F2c_layer_woven, F12_layer_woven, &
                strain, FOS, 2)
    F1c_spar = maxval(FOS)

    ! Tranverse Tensile strength
    stress = [0.0, 1.0, 0.0]
    strain = matmul(Q_lam_woven, stress)
    call TWstrength(E1_layer_woven, E2_layer_woven, theta_layer_woven, G12_layer_woven, nju12_layer_woven, &
                F1t_layer_woven, F2t_layer_woven, F1c_layer_woven, F2c_layer_woven, F12_layer_woven, &
                strain, FOS, 2)
    F2t_spar = maxval(FOS)

    ! Tranverse Compressive strength
    stress = [0.0, -1.0, 0.0]
    strain = matmul(Q_lam_woven, stress)
    call TWstrength(E1_layer_woven, E2_layer_woven, theta_layer_woven, G12_layer_woven, nju12_layer_woven, &
                F1t_layer_woven, F2t_layer_woven, F1c_layer_woven, F2c_layer_woven, F12_layer_woven, &
                strain, FOS, 2)
    F2c_spar = maxval(FOS)

    ! In-plane Shear strength
    stress = [0.0, 0.0, 1.0]
    strain = matmul(Q_lam_woven, stress)
    call TWstrength(E1_layer_woven, E2_layer_woven, theta_layer_woven, G12_layer_woven, nju12_layer_woven, &
                F1t_layer_woven, F2t_layer_woven, F1c_layer_woven, F2c_layer_woven, F12_layer_woven, &
                strain, FOS, 2)
    F12_spar = maxval(FOS)
    
    ! Define properties for spar layers
    E1_layer_spar = E1_woven
    E2_layer_spar = E2_woven
    G12_layer_spar = G12_woven
    nju12_layer_spar = nju12_woven
    rho_spar = rho_282
    t_layer_spar = t_woven
    ! Theta_layer_spar
    do i = 1, n_layer_spar
        theta_layer_spar(i) = (-1)**i*pi/4.0d0
    end do
    t_spar = sum(t_layer_spar)

    ! Module cua spar
    z0 = sum(t_layer_spar) / 2.0
    call laminate(A_lam_spar, D_lam_spar, E1_layer_spar, E2_layer_spar, theta_layer_spar, &
                  t_layer_spar, G12_layer_spar, nju12_layer_spar, z0, size(theta_layer_spar))
    
    Q_lam_spar = A_lam_spar / sum(t_layer_spar)
    S_lam_spar = matinv3(Q_lam_spar)

    Ex_spar = 1.0 / S_lam_spar(1, 1)  ! Gibson Eq. 2.25

    style_ud = 1
    call ply(E1_ud, E2_ud, G12_ud, nju12_ud, sL_tens_ud, sL_compr_ud, sT_tens_ud, &
             sT_compr_ud, sLT_ud, rho_ud, E1f_carb, E2f_carb, G12f_carb, Em, G12m, &
             vf_carb, sLf_tens_carb, sLm_tens, sLm_compr, sLTm, &
             rho_f_carb, rho_m, nju12f_carb, nju12_m, style_ud)

    m_ud = 0.404d0             ! mass of fiber per m^2 in each ply
    Vf_ud = m_ud/rho_f_carb    ! volume of fiber per m^2
    V_ud = Vf_ud/vf_carb       ! total volume of the ply per m^2
    t_ud = V_ud                ! ply thickness
    ! Define properties for cap layers
    E1_layer_cap = E1_ud
    E2_layer_cap = E2_ud
    G12_layer_cap = G12_ud
    theta_layer_cap = 0.0d0
    nju12_layer_cap = nju12_ud
    rho_cap = rho_ud
    t_layer_cap = t_ud

    ! Composite properties (Gibson Chap. 7)
    t_cap = sum(t_layer_cap)

    ! Spar cap calculations
    z0 = h_spar / 2.0 + t_cap
    call laminate(A_lam_cap, D_lam_cap, E1_layer_cap, E2_layer_cap, theta_layer_cap, &
                  t_layer_cap, G12_layer_cap, nju12_layer_cap, z0, size(theta_layer_cap))
    
    Q_lam_cap = A_lam_cap / sum(t_layer_cap)

    ! Tinh do cung cua dam
    w_cap = 55.08d-3
    I_spar = t_spar * h_spar**3 / 12.0  ! Area moment of inertia
    kapa_unit = matmul(matinv3(D_lam_cap),(/1.0d0, 0.0d0, 0.0d0/))  ! Gibson Eq. 7.41
    EI_cap = 1.0 / kapa_unit(1) * w_cap
    EI = 2*Ex_spar * I_spar + 2.0d0*EI_cap  ! Do cung cua dam - Bending rigidity
    ! Tai ap suat, tinh chuyen vi
    w = p * w_cap
    d_max = w * L**4 / 8.0d0 / EI
    ! Output the maximum deflection
    print *, "Maximum deflection (distributed load): ", d_max
    ! Tinh ung suat
    do i = 1, 101
        x(i) = (i-1)*L/100
        M(i) = w*(L-x(i))**2/2.0d0
        kapa(i) = M(i)/EI
        sigma_cap(i) = E1_ud*kapa(i)*(t_cap + h_spar/2)
    end do
    open(10, file = "stress_beamkocore.txt", status = "replace", action = "write")
    do i = 1, 101
        write(10, '(F10.2, 1X, F12.2)') x(i), sigma_cap(i)
    end do
    close(10)
    call system ("python write_stress_result.py")
end program beam_without_core
  