
module LCE_Model

    use iso_fortran_env, only: real64
    implicit none

    private

    integer, parameter, public :: dp = real64

    public :: calculate_r
    public :: build_step_tensor
    public :: matrix_multiply
    public :: print_matrix
    public :: calculate_kinematics
    public :: principal_stretches
    public :: calculate_energy

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
    ! Multiply two 3x3 matrices
    !=========================================================
    subroutine matrix_multiply(A, B, C)

        real(dp), intent(in) :: A(3,3), B(3,3)
        real(dp), intent(out) :: C(3,3)

        integer :: i, j, k

        C = 0.0_dp

        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                end do
            end do
        end do

    end subroutine matrix_multiply


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

        call matrix_multiply(F, Fstar, Ftilde)

        call matrix_multiply(F, invFT, FME)

    end subroutine calculate_kinematics


    !=========================================================
    ! Calculate principal stretches using Jacobi rotations
    !
    ! C = F^T F
    ! lambda_i = sqrt(eigenvalue_i(C))
    !
    ! Stretches are sorted from smallest to largest.
    !=========================================================
    subroutine principal_stretches(F, lambda)

        real(dp), intent(in) :: F(3,3)
        real(dp), intent(out) :: lambda(3)

        real(dp) :: A(3,3)
        real(dp) :: eigenvalues(3)

        real(dp) :: app, aqq, apq
        real(dp) :: akp, akq
        real(dp) :: tau, t, c, s
        real(dp) :: maxoff, temp

        integer :: i, j, k, p, q, iter

        ! Right Cauchy-Green tensor
        A = matmul(transpose(F), F)

        ! Jacobi diagonalization
        do iter = 1, 100

            maxoff = 0.0_dp
            p = 1
            q = 2

            do i = 1, 2
                do j = i+1, 3

                    if (abs(A(i,j)) > maxoff) then
                        maxoff = abs(A(i,j))
                        p = i
                        q = j
                    end if

                end do
            end do

            ! Check convergence
            if (maxoff <= 1.0e-12_dp * &
                max(1.0_dp,maxval(abs(A)))) exit

            app = A(p,p)
            aqq = A(q,q)
            apq = A(p,q)

            tau = (aqq-app)/(2.0_dp*apq)

            t = sign(1.0_dp,tau) / &
                (abs(tau)+sqrt(1.0_dp+tau**2))

            c = 1.0_dp/sqrt(1.0_dp+t**2)
            s = t*c

            ! Update diagonal terms
            A(p,p) = app-t*apq
            A(q,q) = aqq+t*apq

            A(p,q) = 0.0_dp
            A(q,p) = 0.0_dp

            ! Update remaining matrix entries
            do k = 1, 3

                if (k /= p .and. k /= q) then

                    akp = A(k,p)
                    akq = A(k,q)

                    A(k,p) = c*akp-s*akq
                    A(p,k) = A(k,p)

                    A(k,q) = s*akp+c*akq
                    A(q,k) = A(k,q)

                end if

            end do

        end do

        if (iter > 100) then
            error stop "ERROR: Jacobi algorithm did not converge."
        end if

        ! Extract eigenvalues
        do i = 1, 3

            eigenvalues(i) = A(i,i)

            if (eigenvalues(i) < -1.0e-10_dp) then
                error stop "ERROR: Negative eigenvalue in C."
            end if

            eigenvalues(i) = max(eigenvalues(i),0.0_dp)

        end do

        ! Sort eigenvalues in ascending order
        do i = 1, 2
            do j = i+1, 3

                if (eigenvalues(i) > eigenvalues(j)) then

                    temp = eigenvalues(i)
                    eigenvalues(i) = eigenvalues(j)
                    eigenvalues(j) = temp

                end if

            end do
        end do

        ! Principal stretches
        lambda = sqrt(eigenvalues)

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

        character(len=1), intent(out) :: phase

        if (r <= 0.0_dp) error stop "ERROR: r must be positive."
        if (mu < 0.0_dp) error stop "ERROR: mu cannot be negative."
        if (mureg < 0.0_dp) error stop "ERROR: mureg cannot be negative."

        ! Principal stretches of Ftilde
        call principal_stretches(Ftilde, lambda)

        ! Principal stretches of FME
        call principal_stretches(FME, lambdaME)

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

end module LCE_Model
