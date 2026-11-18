
program LCE_Kinematics

    use LCE_Model
    implicit none

    !=========================================================
    ! Input variables
    !=========================================================
    real(dp) :: T0, Tni, T
    real(dp) :: r0, r, rT

    real(dp) :: n0(3)
    real(dp) :: nmag

    real(dp) :: F(3,3)

    !=========================================================
    ! Deformation tensors
    !=========================================================
    real(dp) :: Fstar(3,3)
    real(dp) :: FT(3,3)
    real(dp) :: invFT(3,3)
    real(dp) :: Ftilde(3,3)
    real(dp) :: FME(3,3)

    !=========================================================
    ! Material parameters
    !=========================================================
    real(dp) :: mu
    real(dp) :: mureg

    !=========================================================
    ! Principal stretches
    !=========================================================
    real(dp) :: lambda(3)
    real(dp) :: lambdaME(3)

    !=========================================================
    ! Energy variables
    !=========================================================
    real(dp) :: psi_qc
    real(dp) :: psi_reg
    real(dp) :: psi_total

    character(len=1) :: phase

    integer :: i

    !=========================================================
    ! Appendix Variables
    !=========================================================
    real(dp) :: J
    real(dp) :: Fbar(3,3)
    real(dp) :: FtildeBar(3,3)
    real(dp) :: FMEBar(3,3)

    real(dp) :: Jtilde, JME

    real(dp) :: directions(3,3)
    real(dp) :: directionsME(3,3)

    real(dp) :: lambdaBar(3)
    real(dp) :: lambdaMEBar(3)

    real(dp) :: dWqc(3)
    real(dp) :: dWreg(3)

    real(dp) :: kv
    real(dp) :: sigma(3,3)
    !=========================================================
    ! Read temperature inputs
    !=========================================================
    print *, "=========================================="
    print *, "LIQUID CRYSTAL ELASTOMER MODEL"
    print *, "=========================================="

    print *
    print *, "Enter initial temperature T0:"
    read *, T0

    print *, "Enter nematic-isotropic transition temperature Tni:"
    read *, Tni

    print *, "Enter current temperature T:"
    read *, T

    print *, "Enter initial anisotropy r0:"
    read *, r0

    ! Calculate current anisotropy
    call calculate_r(T, T0, Tni, r0, r)

    rT = r/r0

    !=========================================================
    ! Read initial director
    !=========================================================
    print *
    print *, "Enter initial director n0 (three components):"
    read *, n0(1), n0(2), n0(3)

    ! Normalize director
    nmag = sqrt(sum(n0**2))

    if (nmag <= tiny(1.0_dp)) then
        error stop "ERROR: Director cannot be a zero vector."
    end if

    n0 = n0/nmag

    !=========================================================
    ! Read deformation gradient
    !=========================================================
    print *
    print *, "Enter deformation gradient F row by row:"
    print *, "Enter three values for each row."

    do i = 1, 3
        read *, F(i,1), F(i,2), F(i,3)
    end do

    !=========================================================
    ! Read material parameters
    !=========================================================
    print *
    print *, "Enter shear modulus mu:"
    read *, mu

    print *, "Enter regulation shear modulus mureg:"
    read *, mureg

    !=========================================================
    ! Calculate deformation tensors
    !=========================================================
    call calculate_kinematics(F, n0, r0, r, &
                              Fstar, FT, invFT, Ftilde, FME)

    !=========================================================
    ! Calculate free energies and phases
    !=========================================================
    call calculate_energy(Ftilde, FME, r, mu, mureg, &
                          lambda, lambdaME, phase, &
                          psi_qc, psi_reg, psi_total)

    !=========================================================
! Isochoric decomposition
!=========================================================

call isochoric_decomposition(F, J, Fbar)

call isochoric_decomposition(Ftilde, Jtilde, FtildeBar)

call isochoric_decomposition(FME, JME, FMEBar)

!=========================================================
! Isochoric principal stretches and directions
!=========================================================

call principal_stretches(FtildeBar, lambdaBar, directions)

call principal_stretches(FMEBar, lambdaMEBar, directionsME)

!=========================================================
! First derivatives of free energy
!=========================================================

call calculate_energy_derivatives( &
    lambdaBar, lambdaMEBar, r, mu, mureg, &
    phase, dWqc, dWreg)

!=========================================================
! Bulk modulus
!=========================================================

kv = 500.0_dp*mu

!=========================================================
! Cauchy stress
!=========================================================

call calculate_cauchy_stress( &
    J, lambdaBar, lambdaMEBar, &
    directions, directionsME, &
    dWqc, dWreg, kv, sigma)

    !=========================================================
    ! Display input information
    !=========================================================
    print *
    print *, "=========================================="
    print *, "INPUT PARAMETERS"
    print *, "=========================================="

    print *
    print *, "Initial temperature T0:"
    print *, T0

    print *, "Transition temperature Tni:"
    print *, Tni

    print *, "Current temperature T:"
    print *, T

    print *
    print *, "Initial anisotropy r0:"
    print *, r0

    print *, "Current anisotropy r(T):"
    print *, r

    print *, "Relative anisotropy rT:"
    print *, rT

    print *
    print *, "Normalized director n0:"
    print '(3F14.6)', n0

    !=========================================================
    ! Display deformation tensors
    !=========================================================
    print *
    print *, "=========================================="
    print *, "DEFORMATION TENSORS"
    print *, "=========================================="

    print *
    print *, "Initial spontaneous deformation Fstar:"
    call print_matrix(Fstar)

    print *
    print *, "Thermal deformation gradient FT:"
    call print_matrix(FT)

    print *
    print *, "Inverse thermal deformation gradient:"
    call print_matrix(invFT)

    print *
    print *, "Transformed deformation gradient Ftilde:"
    call print_matrix(Ftilde)

    print *
    print *, "Mechanically transformed deformation FME:"
    call print_matrix(FME)

    !=========================================================
    ! Display principal stretches
    !=========================================================
    print *
    print *, "=========================================="
    print *, "PRINCIPAL STRETCHES"
    print *, "=========================================="

    print *
    print *, "Principal stretches of Ftilde:"
    print '(3F14.6)', lambda

    print *
    print *, "Principal stretches of FME:"
    print '(3F14.6)', lambdaME

    !=========================================================
    ! Display phase and energy results
    !=========================================================
    print *
    print *, "=========================================="
    print *, "LCE FREE ENERGY RESULTS"
    print *, "=========================================="

    print *
    print *, "Shear modulus mu:"
    print *, mu

    print *, "Regulation shear modulus mureg:"
    print *, mureg

    print *
    print *, "Predicted LCE phase: ", phase

    select case (phase)

        case ('L')
            print *, "Phase L: Liquid-like behavior"

        case ('S')
            print *, "Phase S: Solid-like behavior"

        case ('M')
            print *, "Phase M: Mixed liquid/solid behavior"

    end select


    print *
    print *, "=========================================="
    print *, "CAUCHY STRESS TENSOR"
    print *, "=========================================="

    call print_matrix(sigma)

    print *
    print *, "Principal directions of Ftilde:"

    call print_matrix(directions)

    print *
    print *, "Principal directions of FME:"

    call print_matrix(directionsME)

    print *
    print *, "Isochoric principal stretches:"
    print '(3F14.6)', lambdaBar

    print *
    print *, "Isochoric principal stretches of FME:"
    print '(3F14.6)', lambdaMEBar



    print *
    print *, "Quasi-convexified energy:"
    print *, psi_qc

    print *
    print *, "Regulation energy:"
    print *, psi_reg

    print *
    print *, "Total thermo-mechanical energy:"
    print *, psi_total

    print *
    print *, "=========================================="
    print *, "CALCULATION COMPLETE"
    print *, "=========================================="

end program LCE_Kinematics
