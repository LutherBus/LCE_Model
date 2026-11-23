
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
    ! Second derivatives of free energy (Hessians)
    !=========================================================
    real(dp) :: Hqc(3,3)
    real(dp) :: Hreg(3,3)

    !=========================================================
    ! Gamma coefficients for spatial elasticity tensor
    !=========================================================
    real(dp) :: gammaQC(3,3)
    real(dp) :: gammaReg(3,3)

    
    !=========================================================
    ! Spatial elasticity tensor
    !=========================================================
    real(dp) :: CprincipalQC(3,3,3,3)
    real(dp) :: CprincipalReg(3,3,3,3)

    real(dp) :: CmixedQC(3,3,3,3)
    real(dp) :: CmixedReg(3,3,3,3)

    real(dp) :: CisoQC(3,3,3,3)
    real(dp) :: CisoReg(3,3,3,3)

    real(dp) :: Cvol(3,3,3,3)
    real(dp) :: Ctotal(3,3,3,3)


    !=========================================================
    ! Jaumann tangent and Abaqus stiffness matrix
    !=========================================================

    real(dp) :: Cjaumann(3,3,3,3)
    real(dp) :: DDSDDE(6,6)

    !=========================================================
    ! Jaumann tangent verification variables
    !=========================================================

    real(dp) :: Dtest(3,3)
    real(dp) :: hFD

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
! Calculate second derivatives of free energy
!=========================================================
call calculate_second_derivatives( &
    lambdaBar, r, mu, mureg, phase, Hqc, Hreg)
    

!=========================================================
! Calculate gamma coefficients (Equation A.7)
!=========================================================

call calculate_gamma(lambdaBar, dWqc, Hqc, gammaQC)

call calculate_gamma(lambdaMEBar, dWreg, Hreg, gammaReg)

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


    !=========================================================
    ! Display Hessian matrices
    !=========================================================
    print *
    print *, "=========================================="
    print *, "SECOND DERIVATIVES OF FREE ENERGY"
    print *, "=========================================="

    print *
    print *, "Quasi-convexified energy Hessian:"
    call print_matrix(Hqc)

    print *
    print *, "Regulation energy Hessian:"
    call print_matrix(Hreg)


    !=========================================================
    ! Display gamma coefficients
    !=========================================================

    print *
    print *, "=========================================="
    print *, "GAMMA COEFFICIENTS - EQUATION A.7"
    print *, "=========================================="

    print *
    print *, "Quasi-convexified gamma matrix:"
    call print_matrix(gammaQC)

    print *
    print *, "Regulation gamma matrix:"
    call print_matrix(gammaReg)

    print *
    print *, "Quasi-convexified energy:"
    print *, psi_qc

    print *
    print *, "Regulation energy:"
    print *, psi_reg

    print *
    print *, "Total thermo-mechanical energy:"
    print *, psi_total

    call calculate_principal_elasticity( &
        lambdaBar, dWqc, gammaQC, directions, CprincipalQC)

    call calculate_principal_elasticity( &
        lambdaMEBar, dWreg, gammaReg, directionsME, CprincipalReg)

    !=========================================================
    ! Calculate mixed elasticity contributions
    !=========================================================

    call calculate_mixed_elasticity( &
        lambdaBar, dWqc, gammaQC, directions, CmixedQC)

    call calculate_mixed_elasticity( &
        lambdaMEBar, dWreg, gammaReg, directionsME, CmixedReg)

    !=========================================================
    ! Calculate volumetric elasticity
    !=========================================================

    call calculate_volumetric_elasticity(J, kv, Cvol)

    !=========================================================
    ! Assemble complete spatial elasticity tensor
    !=========================================================

    
    !=========================================================
    ! Assemble complete spatial elasticity tensor
    !=========================================================

    CisoQC = CprincipalQC + CmixedQC

    CisoReg = CprincipalReg + CmixedReg

    call calculate_volumetric_elasticity(J, kv, Cvol)

    Ctotal = CisoQC + CisoReg + Cvol

    !=========================================================
    ! Verify tensor assembly
    !=========================================================

    print *
    print *, "COMPLETE SPATIAL ELASTICITY TEST"

    print *, "QC C1212:", CisoQC(1,2,1,2)
    print *, "Reg C1212:", CisoReg(1,2,1,2)
    print *, "Vol C1212:", Cvol(1,2,1,2)

    print *, "Total C1212:", Ctotal(1,2,1,2)

    
    print *
    print *, "PRINCIPAL ELASTICITY TENSOR TEST"

    print *, "QC C1111:", CprincipalQC(1,1,1,1)
    print *, "Reg C1111:", CprincipalReg(1,1,1,1)
    print *, "Reg C1122:", CprincipalReg(1,1,2,2)
    print *, "Reg C1212:", CprincipalReg(1,2,1,2)

    print *
    print *, "=========================================="
    print *, "MIXED ELASTICITY TENSOR TEST"
    print *, "=========================================="

    print *, "QC C1212:", CmixedQC(1,2,1,2)
    print *, "Reg C1212:", CmixedReg(1,2,1,2)

    print *, "QC C1313:", CmixedQC(1,3,1,3)
    print *, "Reg C1313:", CmixedReg(1,3,1,3)

    print *, "Reg C1221:", CmixedReg(1,2,2,1)
    print *, "Reg C1111:", CmixedReg(1,1,1,1)

    print *
    print *, "Mixed tensor symmetry checks:"

    print *, "Minor symmetry error:"
    print *, abs(CmixedReg(1,2,1,2) - CmixedReg(2,1,1,2))

    print *, "Major symmetry error:"
    print *, abs(CmixedReg(1,2,1,3) - CmixedReg(1,3,1,2))

    !=========================================================
    ! Verify combined elasticity tensors
    !=========================================================

    print *
    print *, "=========================================="
    print *, "COMBINED ISOCHORIC ELASTICITY"
    print *, "=========================================="

    print *, "QC total iso C1212:", CisoQC(1,2,1,2)
    print *, "Reg total iso C1212:", CisoReg(1,2,1,2)


    !=========================================================
    ! Volumetric elasticity verification
    !=========================================================

    print *
    print *, "=========================================="
    print *, "VOLUMETRIC ELASTICITY TEST"
    print *, "=========================================="

    print *, "J:", J
    print *, "Cvol 1111:", Cvol(1,1,1,1)
    print *, "Cvol 1122:", Cvol(1,1,2,2)
    print *, "Cvol 1212:", Cvol(1,2,1,2)
    print *, "Cvol 1221:", Cvol(1,2,2,1)

    print *
    print *, "=========================================="
    print *, "COMPLETE SPATIAL ELASTICITY TEST"
    print *, "=========================================="

    print *, "Ctotal 1111:", Ctotal(1,1,1,1)
    print *, "Ctotal 1122:", Ctotal(1,1,2,2)
    print *, "Ctotal 1212:", Ctotal(1,2,1,2)
    
    call check_elasticity_tensor(Ctotal, "Complete elasticity")

    !=========================================================
    ! Calculate Jaumann tangent - Equation A.8
    !=========================================================

    call calculate_jaumann_tangent( &
        Ctotal, sigma, J, Cjaumann)

    !=========================================================
    ! Convert fourth-order tensor to Abaqus Voigt notation
    !=========================================================

    call tensor_to_voigt(Cjaumann, DDSDDE)


    !=========================================================
    ! Finite-difference verification of Jaumann tangent
    !=========================================================

    hFD = 1.0e-5_dp

    ! Symmetric shear perturbation
    Dtest = 0.0_dp
    Dtest(3,3) = 1.0_dp

    print *
    print *, "BASE STATE CAUCHY STRESS:"
    call print_matrix(sigma)

    print *
    print *, "BASE STATE J:"
    print *, J

    print *
    print *, "BASE STATE PHASE:"
    print *, phase

    call verify_jaumann_fd( &
        F, T, T0, Tni, r0, n0, mu, mureg, &
        Cjaumann, hFD, Dtest)


    !=========================================================
    ! Display Jaumann tangent
    !=========================================================

    print *
    print *, "=========================================="
    print *, "JAUMANN TANGENT - EQUATION A.8"
    print *, "=========================================="

    print *, "Cjaumann 1111:", Cjaumann(1,1,1,1)
    print *, "Cjaumann 1122:", Cjaumann(1,1,2,2)
    print *, "Cjaumann 1212:", Cjaumann(1,2,1,2)

    call check_elasticity_tensor(Cjaumann, "Jaumann tangent")

    !=========================================================
    ! Display Abaqus DDSDDE
    !=========================================================

    print *
    print *, "=========================================="
    print *, "ABAQUS DDSDDE MATRIX"
    print *, "=========================================="

    do i = 1,6
        print '(6ES20.10)', DDSDDE(i,:)
    end do

    print *
    print *, "DDSDDE symmetry error:"
    print *, maxval(abs(DDSDDE-transpose(DDSDDE)))

    print *
    print *, "=========================================="
    print *, "CALCULATION COMPLETE"
    print *, "=========================================="

end program LCE_Kinematics
