module ply_module
  implicit none
contains

  subroutine ply(E1c, E2c, G12c, nju12_c, sLc_tens, sLc_compr, sTc_tens, sTc_compr, sLTc, rho_c, &
                 E1f, E2f, G12f, Em, G12m, vf, sLf_tens, sLm_tens, sLm_compr, sLTm, &
                 rho_f, rho_m, nju12_f, nju12_m, style)
    ! Inputs
    real*8, intent(in) :: E1f, E2f, G12f, Em, G12m, vf
    real*8, intent(in) :: sLf_tens, sLm_tens, sLm_compr, sLTm
    real*8, intent(in) :: rho_f, rho_m, nju12_f, nju12_m
    integer, intent(in) :: style

    ! Outputs
    real*8, intent(out) :: E1c, E2c, G12c, nju12_c
    real*8, intent(out) :: sLc_tens, sLc_compr, sTc_tens, sTc_compr, sLTc
    real*8, intent(out) :: rho_c

    ! Local variables
    real*8 :: zeta, eta_E, eta_G, eLf_tens, eLm_tens, phi
    real*8 :: sLc_compr_ext, sLc_compr_shear, sLc_compr_misalign
    real*8 :: d_s_ratio, F

    ! Modulus
    E1c = E1f * vf + Em * (1.0d0 - vf) ! Gibson Eq. 3.23
    nju12_c = nju12_f * vf + nju12_m * (1.0d0 - vf) ! Gibson Eq. 3.41

    zeta = 1.0d0
    eta_E = (E2f - Em) / (E2f + zeta * Em)
    E2c = Em * (1.0d0 + zeta * eta_E * vf) / (1.0d0 - eta_E * vf)
    eta_G = (G12f / G12m - 1.0d0) / (G12f / G12m + zeta)
    G12c = G12m * (1.0d0 + zeta * eta_G * vf) / (1.0d0 - eta_G * vf)

    ! Strength
    eLf_tens = sLf_tens / E1f
    eLm_tens = sLm_tens / Em
    if (eLf_tens < eLm_tens) then
        sLc_tens = eLf_tens * E1c
        if (vf < (sLm_tens - eLf_tens * Em) / (sLf_tens - eLf_tens * Em)) then
            sLc_tens = sLm_tens * (1.0 - vf) ! Gibson Eq. 4.23
        end if
    else
        sLc_tens = eLm_tens * E1c
    end if

    phi = 0.034906
    sLc_compr_ext = 2.0 * vf * sqrt(vf * Em * E1f / 3.0 / (1.0 - vf)) ! Gibson Eq. 4.29
    sLc_compr_shear = G12m / (1.0 - vf) ! Daniel 5.22
    sLc_compr_misalign = sLTm / (sLTm / (G12m / (1.0 - vf)) + phi) ! Gibson Eq. 4.30
    sLc_compr = min(sLc_compr_ext, min(sLc_compr_misalign, sLc_compr_shear))

    ! Transverse strength
    if (style == 1) then
        d_s_ratio = sqrt(vf * 4.0 / 3.14159265)
    else
        d_s_ratio = 0.0
    end if

    F = (d_s_ratio * (Em / E2f - 1.0) + 1.0) ** (-1.0) ! Gibson Eq. 4.38
    sTc_tens = E2c * sLm_tens / Em / F ! Gibson Eq. 4.35
    sTc_compr = E2c * sLm_compr / Em / F ! Gibson Eq. 4.35

    ! Shear strength
    F = (d_s_ratio * (G12m / G12f - 1.0) + 1.0) ** (-1.0) ! Gibson Eq. 4.41
    sLTc = G12c * sLTm / G12m / F ! Gibson Eq. 4.35

    ! Density
    rho_c = rho_f * vf + rho_m * (1.0 - vf)

  end subroutine ply

end module ply_module
