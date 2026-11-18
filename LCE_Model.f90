
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

end module LCE_Model
