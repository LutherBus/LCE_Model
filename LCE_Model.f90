
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


end module LCE_Model
