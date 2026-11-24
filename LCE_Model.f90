
module LCE_Model

    use iso_fortran_env, only: real64
    implicit none

    private

    integer, parameter, public :: dp = real64

    public :: calculate_r
    public :: build_step_tensor
    public :: print_matrix
    public :: calculate_kinematics
    public :: principal_stretches
    public :: calculate_energy
    public :: isochoric_decomposition
    public :: calculate_energy_derivatives
    public :: calculate_cauchy_stress
    public :: calculate_second_derivatives
    public :: calculate_gamma
    public :: fourth_order_dyad
    public :: fourth_order_symdyad
    public :: calculate_principal_elasticity
    public :: calculate_mixed_elasticity
    public :: calculate_volumetric_elasticity
    public :: check_elasticity_tensor
    public :: calculate_jaumann_tangent
    public :: tensor_to_voigt
    public :: evaluate_lce_stress
    public :: verify_jaumann_fd
    public :: verify_full_tangent
    public :: sweep_phase_boundaries
    public :: evaluate_lce_response
    public :: build_sweep_deformation
    public :: verify_boundary_limits

contains

    !=========================================================
    ! Calculate temperature-dependent anisotropy r(T)
    !=========================================================
    subroutine calculate_r(T, T0, Tni, r0, r)

        real(dp), intent(in) :: T, T0, Tni, r0
        real(dp), intent(out) :: r

        real(dp) :: omega, Cr

        if (r0 <= 0.0_dp) error stop "ERROR: r0 must be positive."
        if (Tni <= T0) error stop "ERROR: Tni must exceed T0."

        if (T < T0 .or. T > Tni) then
            error stop "ERROR: Temperature outside T0-Tni range."
        end if

        omega = (T - T0)/(Tni - T0)
        Cr = 1.0_dp - 1.0_dp/r0

        r = 1.0_dp/(Cr*omega + 1.0_dp/r0)

    end subroutine calculate_r


    !=========================================================
    ! Construct a step-length deformation tensor
    !
    ! A = r^(1/3) n0⊗n0 + r^(-1/6)(I - n0⊗n0)
    !=========================================================
    subroutine build_step_tensor(r, n, A)

        real(dp), intent(in) :: r, n(3)
        real(dp), intent(out) :: A(3,3)

        real(dp) :: projector(3,3)
        real(dp) :: identity(3,3)

        integer :: i, j

        if (r <= 0.0_dp) error stop "ERROR: Step ratio must be positive."

        identity = 0.0_dp

        do i = 1, 3
            identity(i,i) = 1.0_dp
        end do

        do i = 1, 3
            do j = 1, 3
                projector(i,j) = n(i)*n(j)
            end do
        end do

        A = r**(1.0_dp/3.0_dp)*projector + &
            r**(-1.0_dp/6.0_dp)*(identity-projector)

    end subroutine build_step_tensor


    !=========================================================
    ! Print a 3x3 matrix
    !=========================================================
    subroutine print_matrix(A)

        real(dp), intent(in) :: A(3,3)

        integer :: i

        do i = 1, 3
            print '(3F14.6)', A(i,:)
        end do

    end subroutine print_matrix


    !=========================================================
    ! Calculate deformation tensors
    !
    ! Fstar = initial spontaneous deformation
    ! FT    = thermal deformation
    ! invFT = inverse thermal deformation
    ! Ftilde = F * Fstar
    ! FME    = F * invFT
    !=========================================================
    subroutine calculate_kinematics(F, n0, r0, r, &
                                    Fstar, FT, invFT, Ftilde, FME)

        real(dp), intent(in) :: F(3,3), n0(3), r0, r

        real(dp), intent(out) :: Fstar(3,3)
        real(dp), intent(out) :: FT(3,3)
        real(dp), intent(out) :: invFT(3,3)
        real(dp), intent(out) :: Ftilde(3,3)
        real(dp), intent(out) :: FME(3,3)

        real(dp) :: rT

        if (r0 <= 0.0_dp) error stop "ERROR: r0 must be positive."

        rT = r/r0

        call build_step_tensor(r0, n0, Fstar)

        call build_step_tensor(rT, n0, FT)

        call build_step_tensor(1.0_dp/rT, n0, invFT)

        Ftilde = matmul(F, Fstar)

        FME = matmul(F, invFT)

    end subroutine calculate_kinematics


    ! !=========================================================
    ! ! Calculate principal stretches using Jacobi rotations
    ! !
    ! ! C = F^T F
    ! ! lambda_i = sqrt(eigenvalue_i(C))
    ! !
    ! ! Stretches are sorted from smallest to largest.
    ! !=========================================================
    ! subroutine principal_stretches(F, lambda)

    !     real(dp), intent(in) :: F(3,3)
    !     real(dp), intent(out) :: lambda(3)

    !     real(dp) :: A(3,3)
    !     real(dp) :: eigenvalues(3)

    !     real(dp) :: app, aqq, apq
    !     real(dp) :: akp, akq
    !     real(dp) :: tau, t, c, s
    !     real(dp) :: maxoff, temp

    !     integer :: i, j, k, p, q, iter

    !     ! Right Cauchy-Green tensor
    !     A = matmul(transpose(F), F)

    !     ! Jacobi diagonalization
    !     do iter = 1, 100

    !         maxoff = 0.0_dp
    !         p = 1
    !         q = 2

    !         do i = 1, 2
    !             do j = i+1, 3

    !                 if (abs(A(i,j)) > maxoff) then
    !                     maxoff = abs(A(i,j))
    !                     p = i
    !                     q = j
    !                 end if

    !             end do
    !         end do

    !         ! Check convergence
    !         if (maxoff <= 1.0e-12_dp * &
    !             max(1.0_dp,maxval(abs(A)))) exit

    !         app = A(p,p)
    !         aqq = A(q,q)
    !         apq = A(p,q)

    !         tau = (aqq-app)/(2.0_dp*apq)

    !         t = sign(1.0_dp,tau) / &
    !             (abs(tau)+sqrt(1.0_dp+tau**2))

    !         c = 1.0_dp/sqrt(1.0_dp+t**2)
    !         s = t*c

    !         ! Update diagonal terms
    !         A(p,p) = app-t*apq
    !         A(q,q) = aqq+t*apq

    !         A(p,q) = 0.0_dp
    !         A(q,p) = 0.0_dp

    !         ! Update remaining matrix entries
    !         do k = 1, 3

    !             if (k /= p .and. k /= q) then

    !                 akp = A(k,p)
    !                 akq = A(k,q)

    !                 A(k,p) = c*akp-s*akq
    !                 A(p,k) = A(k,p)

    !                 A(k,q) = s*akp+c*akq
    !                 A(q,k) = A(k,q)

    !             end if

    !         end do

    !     end do

    !     if (iter > 100) then
    !         error stop "ERROR: Jacobi algorithm did not converge."
    !     end if

    !     ! Extract eigenvalues
    !     do i = 1, 3

    !         eigenvalues(i) = A(i,i)

    !         if (eigenvalues(i) < -1.0e-10_dp) then
    !             error stop "ERROR: Negative eigenvalue in C."
    !         end if

    !         eigenvalues(i) = max(eigenvalues(i),0.0_dp)

    !     end do

    !     ! Sort eigenvalues in ascending order
    !     do i = 1, 2
    !         do j = i+1, 3

    !             if (eigenvalues(i) > eigenvalues(j)) then

    !                 temp = eigenvalues(i)
    !                 eigenvalues(i) = eigenvalues(j)
    !                 eigenvalues(j) = temp

    !             end if

    !         end do
    !     end do

    !     ! Principal stretches
    !     lambda = sqrt(eigenvalues)

    ! end subroutine principal_stretches

    subroutine principal_stretches(F, lambda, directions)

    real(dp), intent(in) :: F(3,3)
    real(dp), intent(out) :: lambda(3)
    real(dp), intent(out) :: directions(3,3)

    real(dp) :: A(3,3), s(3)
    real(dp) :: U(3,3), VT(1,1)
    real(dp) :: wkopt(1)
    real(dp), allocatable :: work(:)

    integer :: info, lwork, i

    external :: dgesvd

    A = F

    ! Workspace query
    lwork = -1

    call dgesvd('A','N',3,3,A,3,s,U,3,VT,1, &
                wkopt,lwork,info)

    if (info /= 0) then
        error stop "ERROR: SVD workspace query failed."
    end if

    lwork = max(1,int(wkopt(1)))
    allocate(work(lwork))

    A = F

    ! Calculate singular values and left singular vectors
    call dgesvd('A','N',3,3,A,3,s,U,3,VT,1, &
                work,lwork,info)

    if (info /= 0) then
        error stop "ERROR: SVD failed."
    end if

    ! LAPACK returns descending singular values.
    ! Reverse order for ascending principal stretches.

    do i = 1,3
        lambda(i) = s(4-i)
        directions(:,i) = U(:,4-i)
    end do

    deallocate(work)

end subroutine principal_stretches

    !=========================================================
    ! Calculate quasi-convexified and regulation energies
    !
    ! Phase L: liquid-like
    ! Phase S: solid-like
    ! Phase M: mixed liquid/solid behavior
    !=========================================================
    subroutine calculate_energy(Ftilde, FME, r, mu, mureg, &
                                lambda, lambdaME, phase, &
                                psi_qc, psi_reg, psi_total)

        real(dp), intent(in) :: Ftilde(3,3), FME(3,3)
        real(dp), intent(in) :: r, mu, mureg

        real(dp), intent(out) :: lambda(3), lambdaME(3)
        real(dp), intent(out) :: psi_qc, psi_reg, psi_total

        real(dp) :: directions(3,3)
        real(dp) :: directionsME(3,3)


        ! Isochoric deformation gradients
        real(dp) :: FtildeBar(3,3)
        real(dp) :: FMEBar(3,3)

        ! Determinants
        real(dp) :: Jtilde, JME


        character(len=1), intent(out) :: phase

        if (r <= 0.0_dp) error stop "ERROR: r must be positive."
        if (mu < 0.0_dp) error stop "ERROR: mu cannot be negative."
        if (mureg < 0.0_dp) error stop "ERROR: mureg cannot be negative."

        call isochoric_decomposition(Ftilde, Jtilde, FtildeBar)
        call isochoric_decomposition(FME, JME, FMEBar)
        
        ! Principal stretches of isochoric deformation
        call principal_stretches(FtildeBar, lambda, directions)
        call principal_stretches(FMEBar, lambdaME, directionsME)   

        ! Avoid division by zero
        if (lambda(1) <= tiny(1.0_dp) .or. &
            lambda(2) <= tiny(1.0_dp)) then

            error stop "ERROR: Deformation gradient is singular."

        end if

        ! Initialize quasi-convexified energy
        psi_qc = 0.0_dp

        !-----------------------------------------
        ! Phase L
        !-----------------------------------------
        if (lambda(1) >= r**(-1.0_dp/6.0_dp)*(1.0_dp-1.0e-7_dp)) then

        phase = 'L'
        psi_qc = 0.0_dp

        !-----------------------------------------
        ! Phase S
        !-----------------------------------------
        else if (lambda(3)/lambda(2) > sqrt(r)) then

            phase = 'S'

            psi_qc = (mu/2.0_dp) * ( &
                r**(1.0_dp/3.0_dp)*lambda(1)**2 + &
                r**(1.0_dp/3.0_dp)*lambda(2)**2 + &
                r**(-2.0_dp/3.0_dp)*lambda(3)**2 - &
                3.0_dp )

        !-----------------------------------------
        ! Phase M
        !-----------------------------------------
        else

            phase = 'M'

            psi_qc = (mu/2.0_dp) * ( &
                r**(1.0_dp/3.0_dp)*lambda(1)**2 + &
                2.0_dp*r**(-1.0_dp/6.0_dp)/lambda(1) - &
                3.0_dp )

        end if

        !-----------------------------------------
        ! Regulation energy
        !-----------------------------------------
        psi_reg = (mureg/2.0_dp) * &
                  (sum(lambdaME**2)-3.0_dp)

        !-----------------------------------------
        ! Total thermo-mechanical energy
        !-----------------------------------------
        psi_total = psi_qc + psi_reg

    end subroutine calculate_energy


subroutine isochoric_decomposition(F, J, Fbar)

    real(dp), intent(in) :: F(3,3)
    real(dp), intent(out) :: J, Fbar(3,3)

    ! Determinant of deformation gradient
    J = F(1,1)*(F(2,2)*F(3,3) - F(2,3)*F(3,2)) &
      - F(1,2)*(F(2,1)*F(3,3) - F(2,3)*F(3,1)) &
      + F(1,3)*(F(2,1)*F(3,2) - F(2,2)*F(3,1))

    if (J <= 0.0_dp) then
        error stop "ERROR: Invalid deformation gradient."
    end if

    ! Isochoric deformation gradient
    Fbar = F / J**(1.0_dp/3.0_dp)

end subroutine isochoric_decomposition

subroutine calculate_energy_derivatives(lambda, lambdaME, &
                                        r, mu, mureg, phase, &
                                        dWqc, dWreg)

    real(dp), intent(in) :: lambda(3), lambdaME(3)
    real(dp), intent(in) :: r, mu, mureg

    character(len=1), intent(in) :: phase

    real(dp), intent(out) :: dWqc(3), dWreg(3)

    dWqc = 0.0_dp
    dWreg = 0.0_dp

    select case (phase)

    case ('L')

        ! Phase L: zero quasi-convexified energy
        dWqc = 0.0_dp

    case ('S')

        ! Phase S: solid-like response
        dWqc(1) = mu*r**(1.0_dp/3.0_dp)*lambda(1)

        dWqc(2) = mu*r**(1.0_dp/3.0_dp)*lambda(2)

        dWqc(3) = mu*r**(-2.0_dp/3.0_dp)*lambda(3)

    case ('M')

        ! Phase M: only minimum stretch contributes
        dWqc(1) = mu*( &
            r**(1.0_dp/3.0_dp)*lambda(1) - &
            r**(-1.0_dp/6.0_dp)/lambda(1)**2 )

        dWqc(2) = 0.0_dp
        dWqc(3) = 0.0_dp

    case default

        error stop "ERROR: Unknown LCE phase."

    end select

    ! Regulation energy derivatives
    dWreg = mureg*lambdaME

end subroutine calculate_energy_derivatives

subroutine calculate_cauchy_stress(J, lambda, lambdaME, &
                                   directions, directionsME, &
                                   dWqc, dWreg, kv, sigma)

    real(dp), intent(in) :: J, kv
    real(dp), intent(in) :: lambda(3), lambdaME(3)
    real(dp), intent(in) :: directions(3,3)
    real(dp), intent(in) :: directionsME(3,3)
    real(dp), intent(in) :: dWqc(3), dWreg(3)

    real(dp), intent(out) :: sigma(3,3)

    real(dp) :: betaQC(3), betaReg(3)
    real(dp) :: sumQC, sumReg
    real(dp) :: dUdJ

    integer :: i, k, a

    if (J <= 0.0_dp) then
        error stop "ERROR: J must be positive."
    end if

    sigma = 0.0_dp

    !-----------------------------------------
    ! Volumetric stress
    ! U = kv*(J-1)^2
    !-----------------------------------------

    dUdJ = 2.0_dp*kv*(J-1.0_dp)

    do i = 1,3
        sigma(i,i) = dUdJ
    end do

    !-----------------------------------------
    ! Equation A.5: beta coefficients
    !-----------------------------------------

    sumQC = sum(lambda*dWqc)
    sumReg = sum(lambdaME*dWreg)

    do a = 1,3

        betaQC(a) = lambda(a)*dWqc(a) &
                    - sumQC/3.0_dp

        betaReg(a) = lambdaME(a)*dWreg(a) &
                     - sumReg/3.0_dp

    end do

    !-----------------------------------------
    ! Equation A.4: isochoric stress
    !-----------------------------------------

    do a = 1,3
        do i = 1,3
            do k = 1,3

                sigma(i,k) = sigma(i,k) + &
                    betaQC(a)/J * &
                    directions(i,a)*directions(k,a)

                sigma(i,k) = sigma(i,k) + &
                    betaReg(a)/J * &
                    directionsME(i,a)*directionsME(k,a)

            end do
        end do
    end do

end subroutine calculate_cauchy_stress


subroutine calculate_second_derivatives(lambda, r, mu, mureg, &
                                        phase, Hqc, Hreg)

    real(dp), intent(in) :: lambda(3)
    real(dp), intent(in) :: r, mu, mureg

    character(len=1), intent(in) :: phase

    real(dp), intent(out) :: Hqc(3,3), Hreg(3,3)

    integer :: i

    ! Initialize Hessian matrices
    Hqc = 0.0_dp
    Hreg = 0.0_dp

    !-----------------------------------------
    ! Quasi-convexified energy Hessian
    !-----------------------------------------

    select case (phase)

    case ('L')

        ! Phase L: Wqc = 0
        Hqc = 0.0_dp

    case ('S')

        ! Phase S
        Hqc(1,1) = mu*r**(1.0_dp/3.0_dp)
        Hqc(2,2) = mu*r**(1.0_dp/3.0_dp)
        Hqc(3,3) = mu*r**(-2.0_dp/3.0_dp)

    case ('M')

        ! Phase M
        Hqc(1,1) = mu*( &
            r**(1.0_dp/3.0_dp) + &
            2.0_dp*r**(-1.0_dp/6.0_dp)/lambda(1)**3 )

    case default

        error stop "ERROR: Unknown LCE phase."

    end select

    !-----------------------------------------
    ! Regulation energy Hessian
    !-----------------------------------------

    do i = 1,3
        Hreg(i,i) = mureg
    end do

end subroutine calculate_second_derivatives


subroutine calculate_gamma(lambda, dW, H, gamma)

    real(dp), intent(in) :: lambda(3)
    real(dp), intent(in) :: dW(3)
    real(dp), intent(in) :: H(3,3)

    real(dp), intent(out) :: gamma(3,3)

    real(dp) :: IM(3,3)
    real(dp) :: rowSum(3)
    real(dp) :: colSum(3)
    real(dp) :: totalSum

    integer :: a, b

    !=========================================================
    ! Construct intermediate matrix IM
    !
    ! IM_ab = lambda_a*lambda_b*H_ab
    !        + delta_ab*lambda_a*dW_a
    !=========================================================

    IM = 0.0_dp

    do a = 1,3
        do b = 1,3

            IM(a,b) = lambda(a)*lambda(b)*H(a,b)

            if (a == b) then
                IM(a,b) = IM(a,b) + lambda(a)*dW(a)
            end if

        end do
    end do

    !=========================================================
    ! Calculate row sums, column sums, and total sum
    !=========================================================

    do a = 1,3
        rowSum(a) = sum(IM(a,:))
        colSum(a) = sum(IM(:,a))
    end do

    totalSum = sum(IM)

    !=========================================================
    ! Calculate gamma coefficients (Equation A.7)
    !=========================================================

    do a = 1,3
        do b = 1,3

            gamma(a,b) = IM(a,b) &
                       - rowSum(a)/3.0_dp &
                       - colSum(b)/3.0_dp &
                       + totalSum/9.0_dp

        end do
    end do

end subroutine calculate_gamma


subroutine fourth_order_dyad(A, B, C)

    real(dp), intent(in) :: A(3,3), B(3,3)
    real(dp), intent(out) :: C(3,3,3,3)

    integer :: i, j, k, l

    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3

                    C(i,j,k,l) = A(i,j)*B(k,l)

                end do
            end do
        end do
    end do

end subroutine fourth_order_dyad


subroutine fourth_order_symdyad(A, B, C)

    real(dp), intent(in) :: A(3,3), B(3,3)
    real(dp), intent(out) :: C(3,3,3,3)

    integer :: i, j, k, l

    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3

                    C(i,j,k,l) = 0.5_dp * ( &
                        A(i,k)*B(j,l) + &
                        A(i,l)*B(j,k) )

                end do
            end do
        end do
    end do

end subroutine fourth_order_symdyad



subroutine calculate_principal_elasticity(lambda, dW, gamma, &
                                         directions, Cprincipal)

    real(dp), intent(in) :: lambda(3), dW(3)
    real(dp), intent(in) :: gamma(3,3)
    real(dp), intent(in) :: directions(3,3)

    real(dp), intent(out) :: Cprincipal(3,3,3,3)

    real(dp) :: beta(3)
    real(dp) :: coefficient
    real(dp) :: projectorA(3,3), projectorB(3,3)
    real(dp) :: dyad(3,3,3,3)

    integer :: a, b, i, k

    ! Calculate beta coefficients (Equation A.5)
    do a = 1,3
        beta(a) = lambda(a)*dW(a) &
                - sum(lambda*dW)/3.0_dp
    end do

    Cprincipal = 0.0_dp

    ! Assemble principal-direction contribution
    do a = 1,3
        do b = 1,3

            coefficient = gamma(a,b)

            if (a == b) then
                coefficient = coefficient - 2.0_dp*beta(a)
            end if

            ! Construct principal projectors
            do i = 1,3
                do k = 1,3
                    projectorA(i,k) = directions(i,a)*directions(k,a)
                    projectorB(i,k) = directions(i,b)*directions(k,b)
                end do
            end do

            ! Fourth-order dyadic product
            call fourth_order_dyad(projectorA, projectorB, dyad)

            Cprincipal = Cprincipal + coefficient*dyad

        end do
    end do

end subroutine calculate_principal_elasticity


subroutine calculate_mixed_elasticity(lambda, dW, gamma, &
                                     directions, Cmixed)

    real(dp), intent(in) :: lambda(3), dW(3)
    real(dp), intent(in) :: gamma(3,3)
    real(dp), intent(in) :: directions(3,3)

    real(dp), intent(out) :: Cmixed(3,3,3,3)

    real(dp) :: beta(3)
    real(dp) :: lambda2(3)
    real(dp) :: q, denominator
    real(dp) :: scale, tolerance
    real(dp) :: meanStress

    integer :: a, b, i, j, k, l

    !=========================================================
    ! Calculate beta coefficients (Equation A.5)
    !=========================================================

    meanStress = sum(lambda*dW)/3.0_dp

    do a = 1,3
        beta(a) = lambda(a)*dW(a) - meanStress
    end do

    ! Squared principal stretches
    lambda2 = lambda**2

    ! Relative tolerance for repeated stretches
    tolerance = 1.0e-8_dp

    Cmixed = 0.0_dp

    !=========================================================
    ! Mixed principal-direction contribution (Equation A.6)
    !=========================================================

    do a = 1,3
        do b = 1,3

            if (a == b) cycle

            denominator = lambda2(b) - lambda2(a)

            scale = max(lambda2(a), lambda2(b))

            !-----------------------------------------
            ! Distinct principal stretches
            !-----------------------------------------

            if (abs(denominator) > tolerance*scale) then

                q = (beta(b)*lambda2(a) - &
                     beta(a)*lambda2(b))/denominator

            !-----------------------------------------
            ! Equal or nearly equal stretches
            !-----------------------------------------

            else

                q = 0.25_dp*( &
                    gamma(a,a) + gamma(b,b) - &
                    gamma(a,b) - gamma(b,a)) &
                    - 0.5_dp*(beta(a) + beta(b))

            end if

            !-----------------------------------------
            ! Assemble fourth-order tensor
            !-----------------------------------------

            do i = 1,3
                do j = 1,3
                    do k = 1,3
                        do l = 1,3

                            Cmixed(i,j,k,l) = &
                                Cmixed(i,j,k,l) + q*( &
                                directions(i,a)*directions(j,b)* &
                                directions(k,a)*directions(l,b) + &
                                directions(i,a)*directions(j,b)* &
                                directions(k,b)*directions(l,a))

                        end do
                    end do
                end do
            end do

        end do
    end do

end subroutine calculate_mixed_elasticity


subroutine calculate_volumetric_elasticity(J, kv, Cvol)

    real(dp), intent(in) :: J, kv
    real(dp), intent(out) :: Cvol(3,3,3,3)

    real(dp) :: identity(3,3)
    real(dp) :: Idyad(3,3,3,3)
    real(dp) :: Isym(3,3,3,3)

    real(dp) :: dUdJ, d2UdJ2
    real(dp) :: coeff1, coeff2

    integer :: i

    !=========================================================
    ! Check input
    !=========================================================

    if (J <= 0.0_dp) then
        error stop "ERROR: J must be positive."
    end if

    if (kv <= 0.0_dp) then
        error stop "ERROR: Bulk modulus must be positive."
    end if

    !=========================================================
    ! Construct identity tensor
    !=========================================================

    identity = 0.0_dp

    do i = 1,3
        identity(i,i) = 1.0_dp
    end do

    !=========================================================
    ! Calculate volumetric energy derivatives
    !=========================================================

    dUdJ = 2.0_dp*kv*(J-1.0_dp)

    d2UdJ2 = 2.0_dp*kv

    !=========================================================
    ! Fourth-order tensor products
    !=========================================================

    call fourth_order_dyad(identity, identity, Idyad)

    call fourth_order_symdyad(identity, identity, Isym)

    !=========================================================
    ! Volumetric elasticity (Equation A.6)
    !=========================================================

    coeff1 = J*dUdJ + J**2*d2UdJ2

    coeff2 = -2.0_dp*J*dUdJ

    Cvol = coeff1*Idyad + coeff2*Isym

end subroutine calculate_volumetric_elasticity


subroutine check_elasticity_tensor(C, label)

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    real(dp), intent(in) :: C(3,3,3,3)
    character(len=*), intent(in) :: label

    real(dp) :: err_minor1, err_minor2, err_major
    integer :: i, j, k, l

    ! Initialize errors
    err_minor1 = 0.0_dp
    err_minor2 = 0.0_dp
    err_major  = 0.0_dp

    !=========================================================
    ! Check for NaN or infinity
    !=========================================================

    if (.not. all(ieee_is_finite(C))) then
        print *, trim(label), ": NON-FINITE VALUES FOUND"
        return
    end if

    !=========================================================
    ! Check fourth-order tensor symmetries
    !=========================================================

    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3

                    ! First minor symmetry: Cijkl = Cjikl
                    err_minor1 = max(err_minor1, &
                        abs(C(i,j,k,l)-C(j,i,k,l)))

                    ! Second minor symmetry: Cijkl = Cijlk
                    err_minor2 = max(err_minor2, &
                        abs(C(i,j,k,l)-C(i,j,l,k)))

                    ! Major symmetry: Cijkl = Cklij
                    err_major = max(err_major, &
                        abs(C(i,j,k,l)-C(k,l,i,j)))

                end do
            end do
        end do
    end do

    !=========================================================
    ! Print verification results
    !=========================================================

    print *
    print *, "=========================================="
    print *, trim(label), " - TENSOR VERIFICATION"
    print *, "=========================================="

    print *, "First minor symmetry error:", err_minor1
    print *, "Second minor symmetry error:", err_minor2
    print *, "Major symmetry error:", err_major

    if (max(err_minor1,err_minor2,err_major) < 1.0e-10_dp) then
        print *, "Symmetry checks: PASSED"
    else
        print *, "Symmetry checks: FAILED"
    end if

end subroutine check_elasticity_tensor


subroutine calculate_jaumann_tangent(Ctotal, sigma, J, Cjaumann)

    real(dp), intent(in) :: Ctotal(3,3,3,3)
    real(dp), intent(in) :: sigma(3,3), J

    real(dp), intent(out) :: Cjaumann(3,3,3,3)

    real(dp) :: identity(3,3)
    real(dp) :: sigmaDyadI(3,3,3,3)
    real(dp) :: IdyadSigma(3,3,3,3)

    integer :: i

    if (J <= 0.0_dp) then
        error stop "ERROR: J must be positive."
    end if

    ! Construct identity tensor
    identity = 0.0_dp

    do i = 1,3
        identity(i,i) = 1.0_dp
    end do

    ! Symmetric fourth-order products
    call fourth_order_symdyad(sigma, identity, sigmaDyadI)
    call fourth_order_symdyad(identity, sigma, IdyadSigma)

    ! Equation A.8
    Cjaumann = Ctotal/J + sigmaDyadI + IdyadSigma

end subroutine calculate_jaumann_tangent


subroutine tensor_to_voigt(Cjaumann, DDSDDE)

    real(dp), intent(in) :: Cjaumann(3,3,3,3)
    real(dp), intent(out) :: DDSDDE(6,6)

    integer :: index1(6), index2(6)
    integer :: a, b
    integer :: i, j, k, l

    ! Abaqus 3D Voigt ordering:
    ! 11, 22, 33, 12, 13, 23

    index1 = [1, 2, 3, 1, 1, 2]
    index2 = [1, 2, 3, 2, 3, 3]

    DDSDDE = 0.0_dp

    do a = 1,6
        do b = 1,6

            i = index1(a)
            j = index2(a)

            k = index1(b)
            l = index2(b)

            DDSDDE(a,b) = Cjaumann(i,j,k,l)

        end do
    end do

end subroutine tensor_to_voigt


subroutine evaluate_lce_stress(F, T, T0, Tni, r0, n0, &
                               mu, mureg, sigma, phase)

    real(dp), intent(in) :: F(3,3), T, T0, Tni
    real(dp), intent(in) :: r0, n0(3), mu, mureg

    real(dp), intent(out) :: sigma(3,3)
    character(len=1), intent(out) :: phase

    real(dp) :: r, kv, J, Jtilde, JME
    real(dp) :: Fstar(3,3), FT(3,3), invFT(3,3)
    real(dp) :: Ftilde(3,3), FME(3,3)
    real(dp) :: Fbar(3,3), FtildeBar(3,3)
    real(dp) :: FMEBar(3,3)

    real(dp) :: lambda(3), lambdaME(3)
    real(dp) :: directions(3,3), directionsME(3,3)
    real(dp) :: dWqc(3), dWreg(3)

    ! Temperature-dependent anisotropy
    call calculate_r(T, T0, Tni, r0, r)

    ! Reference transformations
    call calculate_kinematics(F, n0, r0, r, &
        Fstar, FT, invFT, Ftilde, FME)

    ! Isochoric decomposition
    call isochoric_decomposition(F, J, Fbar)
    call isochoric_decomposition(Ftilde, Jtilde, FtildeBar)
    call isochoric_decomposition(FME, JME, FMEBar)

    ! Principal stretches and directions
    call principal_stretches(FtildeBar, lambda, directions)
    call principal_stretches(FMEBar, lambdaME, directionsME)

    ! Determine LCE phase
    if (lambda(1) >= r**(-1.0_dp/6.0_dp)*(1.0_dp-1.0e-7_dp)) then
        phase = 'L'
    else if (lambda(3)/lambda(2) > sqrt(r)) then
        phase = 'S'
    else
        phase = 'M'
    end if

    ! First derivatives of free energy
    call calculate_energy_derivatives( &
        lambda, lambdaME, r, mu, mureg, &
        phase, dWqc, dWreg)

    ! Bulk modulus
    kv = 500.0_dp*mu

    ! Cauchy stress
    call calculate_cauchy_stress( &
        J, lambda, lambdaME, &
        directions, directionsME, &
        dWqc, dWreg, kv, sigma)

end subroutine evaluate_lce_stress



subroutine verify_jaumann_fd(F, T, T0, Tni, r0, n0, &
                             mu, mureg, Cjaumann, h, D)
    real(dp), intent(in) :: F(3,3), T, T0, Tni, r0, n0(3), mu, mureg
    real(dp), intent(in) :: Cjaumann(3,3,3,3), h, D(3,3)
    real(dp) :: Iden(3,3), Fplus(3,3), Fminus(3,3)
    real(dp) :: sigmaBase(3,3), sigmaPlus(3,3), sigmaMinus(3,3)
    real(dp) :: numericalCauchy(3,3), numericalKirchhoff(3,3)
    real(dp) :: analytical(3,3), analyticalCauchy(3,3)
    real(dp) :: Jbase, Jplus, Jminus, dummyFbar(3,3)
    real(dp) :: absError, relError, scale
    character(len=1) :: phaseBase, phasePlus, phaseMinus
    integer :: i, j, k, l

    if (h <= 0.0_dp) error stop "FD step size must be positive"
    if (maxval(abs(D-transpose(D))) > 1.0e-12_dp) &
        error stop "FD verification requires symmetric D"

    Iden = 0.0_dp
    do i = 1,3
        Iden(i,i) = 1.0_dp
    end do

    Fplus  = matmul(Iden + h*D, F)
    Fminus = matmul(Iden - h*D, F)

    call evaluate_lce_stress(F, T, T0, Tni, r0, n0, &
        mu, mureg, sigmaBase, phaseBase)
    call evaluate_lce_stress(Fplus, T, T0, Tni, r0, n0, &
        mu, mureg, sigmaPlus, phasePlus)
    call evaluate_lce_stress(Fminus, T, T0, Tni, r0, n0, &
        mu, mureg, sigmaMinus, phaseMinus)

    call isochoric_decomposition(F, Jbase, dummyFbar)
    call isochoric_decomposition(Fplus, Jplus, dummyFbar)
    call isochoric_decomposition(Fminus, Jminus, dummyFbar)

    numericalCauchy = (sigmaPlus-sigmaMinus)/(2.0_dp*h)
    numericalKirchhoff = (Jplus*sigmaPlus-Jminus*sigmaMinus) &
                         /(2.0_dp*h*Jbase)

    ! Calculate BEFORE any verification output.
    analytical = 0.0_dp
    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3
                    analytical(i,j) = analytical(i,j) + &
                        Cjaumann(i,j,k,l)*D(k,l)
                end do
            end do
        end do
    end do
    analyticalCauchy = analytical - sigmaBase* &
                        (D(1,1)+D(2,2)+D(3,3))

    absError = maxval(abs(numericalKirchhoff-analytical))
    scale = max(maxval(abs(numericalKirchhoff)),maxval(abs(analytical)))
    relError = absError/max(scale,1.0e-12_dp)
    print *
    print *, "========== ABAQUS KIRCHHOFF FD CHECK =========="
    print *, "Step size:", h
    print *, "Phases (base,+,-):", phaseBase, phasePlus, phaseMinus
    print *, "Numerical (delta tau / J):"
    call print_matrix(numericalKirchhoff)
    print *, "Analytical (Cjaumann:D):"
    call print_matrix(analytical)
    print *, "Maximum absolute error:", absError
    print *, "Relative error:", relError

    absError = maxval(abs(numericalCauchy-analyticalCauchy))
    scale = max(maxval(abs(numericalCauchy)),maxval(abs(analyticalCauchy)))
    relError = absError/max(scale,1.0e-12_dp)
    print *
    print *, "========== DIRECT CAUCHY FD CHECK =========="
    print *, "Numerical delta sigma:"
    call print_matrix(numericalCauchy)
    print *, "Analytical (Cjaumann:D - sigma*tr(D)):"
    call print_matrix(analyticalCauchy)
    print *, "Maximum absolute error:", absError
    print *, "Relative error:", relError
end subroutine verify_jaumann_fd



subroutine verify_full_tangent(F, T, T0, Tni, r0, n0, &
                               mu, mureg, DDSDDE, h)

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    real(dp), intent(in) :: F(3,3), T, T0, Tni
    real(dp), intent(in) :: r0, n0(3), mu, mureg
    real(dp), intent(in) :: DDSDDE(6,6)
    real(dp), intent(in) :: h

    real(dp) :: Iden(3,3), D(3,3)
    real(dp) :: Fplus(3,3), Fminus(3,3)
    real(dp) :: FbarBase(3,3)
    real(dp) :: FbarPlus(3,3), FbarMinus(3,3)

    real(dp) :: sigmaBase(3,3)
    real(dp) :: sigmaPlus(3,3), sigmaMinus(3,3)

    real(dp) :: numerical(3,3)
    real(dp) :: DDSDDE_FD(6,6)
    real(dp) :: errorMatrix(6,6)

    real(dp) :: Jbase, Jplus, Jminus
    real(dp) :: absError, relError, scale
    real(dp) :: colError, colScale

    integer :: index1(6), index2(6)
    integer :: i, a, b, p, q

    character(len=1) :: phaseBase
    character(len=1) :: phasePlus, phaseMinus

    logical :: phaseChanged

    if (h <= 0.0_dp) then
        error stop "ERROR: FD step size must be positive."
    end if

    !=========================================================
    ! Abaqus Voigt ordering
    !=========================================================

    index1 = [1, 2, 3, 1, 1, 2]
    index2 = [1, 2, 3, 2, 3, 3]

    !=========================================================
    ! Identity matrix
    !=========================================================

    Iden = 0.0_dp

    do i = 1,3
        Iden(i,i) = 1.0_dp
    end do

    !=========================================================
    ! Base-state determinant and stress
    !=========================================================

    call isochoric_decomposition(F, Jbase, FbarBase)

    call evaluate_lce_stress(F, T, T0, Tni, r0, n0, &
        mu, mureg, sigmaBase, phaseBase)

    DDSDDE_FD = 0.0_dp
    phaseChanged = .false.

    print *
    print *, "=========================================="
    print *, "FULL DDSDDE FINITE-DIFFERENCE VERIFICATION"
    print *, "=========================================="

    print *, "Base phase:", phaseBase
    print *, "Base J:", Jbase
    print *, "Step size:", h

    !=========================================================
    ! Loop over all six independent strain perturbations
    !=========================================================

    do b = 1,6

        D = 0.0_dp

        p = index1(b)
        q = index2(b)

        ! Normal perturbation
        if (p == q) then

            D(p,q) = 1.0_dp

        ! Engineering shear perturbation
        else

            D(p,q) = 0.5_dp
            D(q,p) = 0.5_dp

        end if

        ! Perturbed deformation gradients
        Fplus  = matmul(Iden + h*D, F)
        Fminus = matmul(Iden - h*D, F)

        ! Perturbed determinants
        call isochoric_decomposition( &
            Fplus, Jplus, FbarPlus)

        call isochoric_decomposition( &
            Fminus, Jminus, FbarMinus)

        ! Perturbed Cauchy stresses
        call evaluate_lce_stress( &
            Fplus, T, T0, Tni, r0, n0, &
            mu, mureg, sigmaPlus, phasePlus)

        call evaluate_lce_stress( &
            Fminus, T, T0, Tni, r0, n0, &
            mu, mureg, sigmaMinus, phaseMinus)

        ! Check for phase-boundary crossings
        print *, "Column:", b, &
                 " Phases:", phasePlus, phaseMinus

        if (phasePlus /= phaseBase .or. &
            phaseMinus /= phaseBase) then

            phaseChanged = .true.

            print *, "WARNING: Phase boundary crossed."

        end if

        !=====================================================
        ! Numerical Abaqus tangent:
        ! derivative of Kirchhoff stress divided by base J
        !=====================================================

        numerical = &
            (Jplus*sigmaPlus - Jminus*sigmaMinus) / &
            (2.0_dp*h*Jbase)

        ! Store derivative in column b
        do a = 1,6

            DDSDDE_FD(a,b) = &
                numerical(index1(a),index2(a))

        end do

    end do

    !=========================================================
    ! Check for invalid numerical values
    !=========================================================

    if (.not. all(ieee_is_finite(DDSDDE_FD)) .or. &
        .not. all(ieee_is_finite(DDSDDE))) then

        print *, "FAILED: Non-finite tangent components."
        return

    end if

    !=========================================================
    ! Calculate error matrix
    !=========================================================

    errorMatrix = DDSDDE_FD - DDSDDE

    absError = maxval(abs(errorMatrix))

    scale = max(maxval(abs(DDSDDE_FD)), &
                maxval(abs(DDSDDE)), 1.0e-12_dp)

    relError = absError/scale

    !=========================================================
    ! Print analytical matrix
    !=========================================================

    print *
    print *, "Analytical DDSDDE:"

    do i = 1,6
        print '(6ES18.8)', DDSDDE(i,:)
    end do

    !=========================================================
    ! Print numerical matrix
    !=========================================================

    print *
    print *, "Numerical DDSDDE:"

    do i = 1,6
        print '(6ES18.8)', DDSDDE_FD(i,:)
    end do

    !=========================================================
    ! Print difference matrix
    !=========================================================

    print *
    print *, "Difference (Numerical - Analytical):"

    do i = 1,6
        print '(6ES18.8)', errorMatrix(i,:)
    end do

    !=========================================================
    ! Print errors for each column
    !=========================================================

    print *
    print *, "COLUMN-BY-COLUMN ERRORS"

    do b = 1,6

        colError = maxval(abs(errorMatrix(:,b)))

        colScale = max( &
            maxval(abs(DDSDDE_FD(:,b))), &
            maxval(abs(DDSDDE(:,b))), 1.0e-12_dp)

        print *, "Column:", b
        print *, "Absolute error:", colError
        print *, "Relative error:", colError/colScale

    end do

    !=========================================================
    ! Overall verification results
    !=========================================================

    print *
    print *, "Maximum absolute error:", absError
    print *, "Overall relative error:", relError

    print *, "Numerical symmetry error:", &
        maxval(abs(DDSDDE_FD-transpose(DDSDDE_FD)))

    if (phaseChanged) then

        print *, "WARNING: Phase change detected."
        print *, "Smooth tangent comparison may be invalid."

    end if

end subroutine verify_full_tangent



subroutine sweep_phase_boundaries(T, T0, Tni, r0, n0, &
                                  mu, mureg, sweepType)

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    real(dp), intent(in) :: T, T0, Tni
    real(dp), intent(in) :: r0, n0(3)
    real(dp), intent(in) :: mu, mureg

    integer, intent(in) :: sweepType

    real(dp) :: r, kv, s, sStart, sEnd

    real(dp) :: F(3,3), Ftarget(3,3)
    real(dp) :: Fstar(3,3), invFstar(3,3)
    real(dp) :: FT(3,3), invFT(3,3)
    real(dp) :: Ftilde(3,3), FME(3,3)

    real(dp) :: Fbar(3,3)
    real(dp) :: FtildeBar(3,3), FMEBar(3,3)

    real(dp) :: J, Jtilde, JME

    real(dp) :: lambda(3), lambdaME(3)
    real(dp) :: lambdaBar(3), lambdaMEBar(3)

    real(dp) :: directions(3,3)
    real(dp) :: directionsME(3,3)

    real(dp) :: psi_qc, psi_reg, psi_total
    real(dp) :: psi_vol

    real(dp) :: dWqc(3), dWreg(3)
    real(dp) :: Hqc(3,3), Hreg(3,3)

    real(dp) :: gammaQC(3,3), gammaReg(3,3)

    real(dp) :: CprincipalQC(3,3,3,3)
    real(dp) :: CprincipalReg(3,3,3,3)

    real(dp) :: CmixedQC(3,3,3,3)
    real(dp) :: CmixedReg(3,3,3,3)

    real(dp) :: Cvol(3,3,3,3)
    real(dp) :: Ctotal(3,3,3,3)
    real(dp) :: Cjaumann(3,3,3,3)

    real(dp) :: sigma(3,3)
    real(dp) :: DDSDDE(6,6)

    integer :: step, nSteps, fileUnit

    character(len=1) :: phase, previousPhase
    character(len=80) :: filename

    !=========================================================
    ! Material parameters
    !=========================================================

    call calculate_r(T, T0, Tni, r0, r)

    kv = 500.0_dp*mu

    call build_step_tensor(1.0_dp/r0, n0, invFstar)

    nSteps = 2000

    !=========================================================
    ! Select deformation path
    !=========================================================

    select case (sweepType)

    case (1)

        sStart = 1.0_dp
        sEnd = 0.35_dp

        filename = "phase_sweep_LMS.csv"

    case (2)

        sStart = 1.0_dp
        sEnd = 2.5_dp

        filename = "phase_sweep_LS.csv"

    case default

        error stop "ERROR: Invalid sweep type."

    end select

    !=========================================================
    ! Open output file
    !=========================================================

    open(newunit=fileUnit, file=trim(filename), &
         status="replace", action="write")

    write(fileUnit,'(A)') &
        "s,phase,J,lambda1,lambda2,lambda3," // &
        "Wqc,Wreg,Wvol,Wtotal," // &
        "sigma11,sigma22,sigma33,sigma12," // &
        "C11,C22,C33,C12,C44,C55,C66"

    previousPhase = ' '

    !=========================================================
    ! Deformation sweep
    !=========================================================

    do step = 0,nSteps

        s = sStart + (sEnd-sStart)* &
            real(step,dp)/real(nSteps,dp)

        Ftarget = 0.0_dp

        select case (sweepType)

        case (1)

            Ftarget(1,1) = 1.0_dp/s
            Ftarget(2,2) = s
            Ftarget(3,3) = 1.0_dp

        case (2)

            Ftarget(1,1) = 1.0_dp/sqrt(s)
            Ftarget(2,2) = 1.0_dp/sqrt(s)
            Ftarget(3,3) = s

        end select

        ! Convert transformed deformation into F
        F = matmul(Ftarget, invFstar)

        !=====================================================
        ! Kinematics
        !=====================================================

        call calculate_kinematics(F, n0, r0, r, &
            Fstar, FT, invFT, Ftilde, FME)

        call isochoric_decomposition(F, J, Fbar)

        call isochoric_decomposition( &
            Ftilde, Jtilde, FtildeBar)

        call isochoric_decomposition( &
            FME, JME, FMEBar)

        !=====================================================
        ! Energy and phase
        !=====================================================

        call calculate_energy(Ftilde, FME, r, mu, mureg, &
            lambda, lambdaME, phase, &
            psi_qc, psi_reg, psi_total)

        psi_vol = kv*(J-1.0_dp)**2

        !=====================================================
        ! Principal stretches and directions
        !=====================================================

        call principal_stretches( &
            FtildeBar, lambdaBar, directions)

        call principal_stretches( &
            FMEBar, lambdaMEBar, directionsME)

        !=====================================================
        ! Energy derivatives
        !=====================================================

        call calculate_energy_derivatives( &
            lambdaBar, lambdaMEBar, r, mu, mureg, &
            phase, dWqc, dWreg)

        call calculate_second_derivatives( &
            lambdaBar, r, mu, mureg, phase, Hqc, Hreg)

        !=====================================================
        ! Cauchy stress
        !=====================================================

        call calculate_cauchy_stress( &
            J, lambdaBar, lambdaMEBar, &
            directions, directionsME, &
            dWqc, dWreg, kv, sigma)

        !=====================================================
        ! Gamma coefficients
        !=====================================================

        call calculate_gamma( &
            lambdaBar, dWqc, Hqc, gammaQC)

        call calculate_gamma( &
            lambdaMEBar, dWreg, Hreg, gammaReg)

        !=====================================================
        ! Principal elasticity
        !=====================================================

        call calculate_principal_elasticity( &
            lambdaBar, dWqc, gammaQC, &
            directions, CprincipalQC)

        call calculate_principal_elasticity( &
            lambdaMEBar, dWreg, gammaReg, &
            directionsME, CprincipalReg)

        !=====================================================
        ! Mixed elasticity
        !=====================================================

        call calculate_mixed_elasticity( &
            lambdaBar, dWqc, gammaQC, &
            directions, CmixedQC)

        call calculate_mixed_elasticity( &
            lambdaMEBar, dWreg, gammaReg, &
            directionsME, CmixedReg)

        !=====================================================
        ! Volumetric and total elasticity
        !=====================================================

        call calculate_volumetric_elasticity(J, kv, Cvol)

        Ctotal = CprincipalQC + CmixedQC + &
                 CprincipalReg + CmixedReg + Cvol

        !=====================================================
        ! Jaumann tangent and DDSDDE
        !=====================================================

        call calculate_jaumann_tangent( &
            Ctotal, sigma, J, Cjaumann)

        call tensor_to_voigt(Cjaumann, DDSDDE)

        !=====================================================
        ! Check numerical stability
        !=====================================================

        if (.not. all(ieee_is_finite(sigma)) .or. &
            .not. all(ieee_is_finite(DDSDDE)) .or. &
            .not. ieee_is_finite(psi_total)) then

            print *, "ERROR: Non-finite values at s =", s

            close(fileUnit)

            error stop "Phase sweep failed."

        end if

        !=====================================================
        ! Detect phase transitions
        !=====================================================

        if (phase /= previousPhase) then

            print *, "s =", s, " Phase =", phase

            previousPhase = phase

        end if

        !=====================================================
        ! Write results to CSV
        !=====================================================

        write(fileUnit,'(ES24.16,A,A1,19(A,ES24.16))') &
            s, ",", phase, &
            ",", J, &
            ",", lambdaBar(1), &
            ",", lambdaBar(2), &
            ",", lambdaBar(3), &
            ",", psi_qc, &
            ",", psi_reg, &
            ",", psi_vol, &
            ",", psi_total + psi_vol, &
            ",", sigma(1,1), &
            ",", sigma(2,2), &
            ",", sigma(3,3), &
            ",", sigma(1,2), &
            ",", DDSDDE(1,1), &
            ",", DDSDDE(2,2), &
            ",", DDSDDE(3,3), &
            ",", DDSDDE(1,2), &
            ",", DDSDDE(4,4), &
            ",", DDSDDE(5,5), &
            ",", DDSDDE(6,6)

    end do

    close(fileUnit)

    print *
    print *, "Phase sweep complete."
    print *, "Results saved to: ", trim(filename)

end subroutine sweep_phase_boundaries


subroutine evaluate_lce_response(F, T, T0, Tni, r0, n0, &
                                 mu, mureg, phase, &
                                 energy, sigma, DDSDDE)

    real(dp), intent(in) :: F(3,3)
    real(dp), intent(in) :: T, T0, Tni, r0
    real(dp), intent(in) :: n0(3), mu, mureg

    character(len=1), intent(out) :: phase

    real(dp), intent(out) :: energy
    real(dp), intent(out) :: sigma(3,3)
    real(dp), intent(out) :: DDSDDE(6,6)

    real(dp) :: r, kv
    real(dp) :: J, Jtilde, JME

    real(dp) :: Fstar(3,3), FT(3,3), invFT(3,3)
    real(dp) :: Ftilde(3,3), FME(3,3)

    real(dp) :: Fbar(3,3)
    real(dp) :: FtildeBar(3,3), FMEBar(3,3)

    real(dp) :: lambda(3), lambdaME(3)
    real(dp) :: lambdaBar(3), lambdaMEBar(3)

    real(dp) :: directions(3,3), directionsME(3,3)

    real(dp) :: psi_qc, psi_reg, psi_total

    real(dp) :: dWqc(3), dWreg(3)
    real(dp) :: Hqc(3,3), Hreg(3,3)
    real(dp) :: gammaQC(3,3), gammaReg(3,3)

    real(dp) :: CprincipalQC(3,3,3,3)
    real(dp) :: CprincipalReg(3,3,3,3)
    real(dp) :: CmixedQC(3,3,3,3)
    real(dp) :: CmixedReg(3,3,3,3)

    real(dp) :: Cvol(3,3,3,3)
    real(dp) :: Ctotal(3,3,3,3)
    real(dp) :: Cjaumann(3,3,3,3)

    ! Temperature-dependent material parameters
    call calculate_r(T, T0, Tni, r0, r)

    kv = 500.0_dp*mu

    ! Kinematics
    call calculate_kinematics(F, n0, r0, r, &
        Fstar, FT, invFT, Ftilde, FME)

    ! Isochoric decomposition
    call isochoric_decomposition(F, J, Fbar)

    call isochoric_decomposition( &
        Ftilde, Jtilde, FtildeBar)

    call isochoric_decomposition( &
        FME, JME, FMEBar)

    ! Free energy and phase selection
    call calculate_energy(Ftilde, FME, r, mu, mureg, &
        lambda, lambdaME, phase, &
        psi_qc, psi_reg, psi_total)

    energy = psi_total + kv*(J-1.0_dp)**2

    ! Principal stretches and spatial directions
    call principal_stretches( &
        FtildeBar, lambdaBar, directions)

    call principal_stretches( &
        FMEBar, lambdaMEBar, directionsME)

    ! First and second energy derivatives
    call calculate_energy_derivatives( &
        lambdaBar, lambdaMEBar, r, mu, mureg, &
        phase, dWqc, dWreg)

    call calculate_second_derivatives( &
        lambdaBar, r, mu, mureg, phase, Hqc, Hreg)

    ! Cauchy stress
    call calculate_cauchy_stress( &
        J, lambdaBar, lambdaMEBar, &
        directions, directionsME, &
        dWqc, dWreg, kv, sigma)

    ! Gamma coefficients
    call calculate_gamma( &
        lambdaBar, dWqc, Hqc, gammaQC)

    call calculate_gamma( &
        lambdaMEBar, dWreg, Hreg, gammaReg)

    ! Principal elasticity
    call calculate_principal_elasticity( &
        lambdaBar, dWqc, gammaQC, &
        directions, CprincipalQC)

    call calculate_principal_elasticity( &
        lambdaMEBar, dWreg, gammaReg, &
        directionsME, CprincipalReg)

    ! Mixed elasticity
    call calculate_mixed_elasticity( &
        lambdaBar, dWqc, gammaQC, &
        directions, CmixedQC)

    call calculate_mixed_elasticity( &
        lambdaMEBar, dWreg, gammaReg, &
        directionsME, CmixedReg)

    ! Volumetric elasticity
    call calculate_volumetric_elasticity(J, kv, Cvol)

    ! Complete spatial elasticity
    Ctotal = CprincipalQC + CmixedQC + &
             CprincipalReg + CmixedReg + Cvol

    ! Abaqus tangent
    call calculate_jaumann_tangent( &
        Ctotal, sigma, J, Cjaumann)

    call tensor_to_voigt(Cjaumann, DDSDDE)

end subroutine evaluate_lce_response



subroutine build_sweep_deformation(s, sweepType, &
                                   r0, n0, F)

    real(dp), intent(in) :: s, r0, n0(3)
    integer, intent(in) :: sweepType

    real(dp), intent(out) :: F(3,3)

    real(dp) :: Ftarget(3,3)
    real(dp) :: invFstar(3,3)

    if (s <= 0.0_dp) then
        error stop "ERROR: s must be positive."
    end if

    call build_step_tensor(1.0_dp/r0, n0, invFstar)

    Ftarget = 0.0_dp

    select case (sweepType)

    case (1)

        ! L-M-S deformation path
        Ftarget(1,1) = 1.0_dp/s
        Ftarget(2,2) = s
        Ftarget(3,3) = 1.0_dp

    case (2)

        ! L-S deformation path
        Ftarget(1,1) = 1.0_dp/sqrt(s)
        Ftarget(2,2) = 1.0_dp/sqrt(s)
        Ftarget(3,3) = s

    case default

        error stop "ERROR: Invalid sweep type."

    end select

    F = matmul(Ftarget, invFstar)

end subroutine build_sweep_deformation


subroutine verify_boundary_limits(T, T0, Tni, r0, n0, &
                                   mu, mureg, boundaryType)

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    real(dp), intent(in) :: T, T0, Tni, r0
    real(dp), intent(in) :: n0(3), mu, mureg

    integer, intent(in) :: boundaryType

    real(dp) :: r, sc, offset, s, h

    real(dp) :: F(3,3), F1(3,3), F2(3,3)

    real(dp) :: energyMinus, energyPlus
    real(dp) :: sigmaMinus(3,3), sigmaPlus(3,3)
    real(dp) :: Cminus(6,6), Cplus(6,6)

    real(dp) :: energyBase, energy1, energy2
    real(dp) :: sigmaBase(3,3)
    real(dp) :: sigma1(3,3), sigma2(3,3)

    real(dp) :: Cbase(6,6), C1(6,6), C2(6,6)

    real(dp) :: Dvoigt(6)
    real(dp) :: analytical(6), numerical(6)
    real(dp) :: derivative(3,3)

    real(dp) :: absError, relError, scale

    integer :: sweepType, iteration, side
    integer :: a

    integer :: index1(6), index2(6)

    character(len=1) :: phaseMinus, phasePlus
    character(len=1) :: phaseBase, phase1, phase2

    character(len=20) :: boundaryName

    call calculate_r(T, T0, Tni, r0, r)

    index1 = [1, 2, 3, 1, 1, 2]
    index2 = [1, 2, 3, 2, 3, 3]

    ! Select theoretical phase boundary
    select case (boundaryType)

    case (1)

        sc = r**(-1.0_dp/6.0_dp)
        sweepType = 1
        boundaryName = "L-M"

    case (2)

        sc = r**(-0.5_dp)
        sweepType = 1
        boundaryName = "M-S"

    case (3)

        sc = r**(1.0_dp/3.0_dp)
        sweepType = 2
        boundaryName = "L-S"

    case default

        error stop "ERROR: Invalid boundary type."

    end select

    print *
    print *, "=========================================="
    print *, "BOUNDARY LIMIT TEST: ", trim(boundaryName)
    print *, "=========================================="

    print *, "Critical s:", sc

    ! Approach boundary from both sides
    do iteration = 2,6

        offset = 10.0_dp**(-iteration)

        ! Negative side
        call build_sweep_deformation( &
            sc-offset, sweepType, r0, n0, F)

        call evaluate_lce_response( &
            F, T, T0, Tni, r0, n0, mu, mureg, &
            phaseMinus, energyMinus, sigmaMinus, Cminus)

        ! Positive side
        call build_sweep_deformation( &
            sc+offset, sweepType, r0, n0, F)

        call evaluate_lce_response( &
            F, T, T0, Tni, r0, n0, mu, mureg, &
            phasePlus, energyPlus, sigmaPlus, Cplus)

        ! Continuity checks
        print *
        print *, "Offset:", offset

        print *, "Phases (-,+):", phaseMinus, phasePlus

        print *, "Energy difference:", &
            abs(energyPlus-energyMinus)

        print *, "Maximum stress difference:", &
            maxval(abs(sigmaPlus-sigmaMinus))

        print *, "Maximum tangent difference:", &
            maxval(abs(Cplus-Cminus))

        if (.not. all(ieee_is_finite(Cminus)) .or. &
            .not. all(ieee_is_finite(Cplus))) then

            error stop "ERROR: Non-finite boundary tangent."

        end if

        !=================================================
        ! One-sided tangent verification on both sides
        !=================================================

        do side = -1,1,2

            s = sc + real(side,dp)*offset

            ! Small step directed away from boundary
            h = 0.1_dp*offset

            call build_sweep_deformation( &
                s, sweepType, r0, n0, F)

            call build_sweep_deformation( &
                s+real(side,dp)*h, &
                sweepType, r0, n0, F1)

            call build_sweep_deformation( &
                s+2.0_dp*real(side,dp)*h, &
                sweepType, r0, n0, F2)

            call evaluate_lce_response( &
                F, T, T0, Tni, r0, n0, mu, mureg, &
                phaseBase, energyBase, sigmaBase, Cbase)

            call evaluate_lce_response( &
                F1, T, T0, Tni, r0, n0, mu, mureg, &
                phase1, energy1, sigma1, C1)

            call evaluate_lce_response( &
                F2, T, T0, Tni, r0, n0, mu, mureg, &
                phase2, energy2, sigma2, C2)

            ! Ensure all states belong to same phase
            if (phaseBase /= phase1 .or. &
                phaseBase /= phase2) then

                print *, "WARNING: Phase crossing."
                cycle

            end if

            ! Second-order one-sided derivative
            derivative = real(side,dp)* &
                (-3.0_dp*sigmaBase + &
                 4.0_dp*sigma1 - sigma2) / &
                (2.0_dp*h)

            ! Spatial deformation rate along sweep path
            Dvoigt = 0.0_dp

            select case (sweepType)

            case (1)

                Dvoigt(1) = -1.0_dp/s
                Dvoigt(2) =  1.0_dp/s

            case (2)

                Dvoigt(1) = -0.5_dp/s
                Dvoigt(2) = -0.5_dp/s
                Dvoigt(3) =  1.0_dp/s

            end select

            ! Analytical stress derivative
            analytical = matmul(Cbase, Dvoigt)

            ! Numerical stress derivative in Voigt form
            do a = 1,6

                numerical(a) = &
                    derivative(index1(a),index2(a))

            end do

            ! Calculate error
            absError = maxval(abs(numerical-analytical))

            scale = max(maxval(abs(numerical)), &
                        maxval(abs(analytical)), &
                        1.0e-12_dp)

            relError = absError/scale

            print *, "Side:", side, " Phase:", phaseBase

            print *, "FD absolute error:", absError
            print *, "FD relative error:", relError

        end do

    end do

end subroutine verify_boundary_limits



end module LCE_Model
