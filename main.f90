program LCE_main
    use LCE_Model, only: dp, evaluate_lce_response
    implicit none
    real(dp) :: T, T0, Tni, r0, n0(3), mu, mureg, F(3,3)
    real(dp) :: energy, sigma(3,3), DDSDDE(6,6)
    character(len=1) :: phase
    integer :: i
    print *, "Enter T0, Tni, T (one per line):"
    read *, T0
    read *, Tni
    read *, T
    print *, "Enter r0:"
    read *, r0
    print *, "Enter n0 (three components):"
    read *, n0
    print *, "Enter F row by row:"
    do i=1,3
        read *, F(i,:)
    end do
    print *, "Enter mu and mureg (one per line):"
    read *, mu
    read *, mureg
    call evaluate_lce_response(F, T, T0, Tni, r0, n0, &
                               mu, mureg, phase, energy, sigma, DDSDDE)
    print *, "Phase: ", phase
    print *, "Total free energy: ", energy
    print *, "Cauchy stress:"
    do i=1,3
        print '(3ES20.10)', sigma(i,:)
    end do
    print *, "DDSDDE:"
    do i=1,6
        print '(6ES20.10)', DDSDDE(i,:)
    end do
end program LCE_main
