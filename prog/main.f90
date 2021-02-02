program main
implicit none
    complex :: H2(2, 2), H3(3, 3)
    integer :: i

    call H_2lvl(H2, 1., 1., 1.)
    call H_3lvl(H3, 1., 1., 1., 1., 1., 1.)

    do i=1, size(H3, 1)
        write(*, *) H3(i, :)
    enddo

end program main


subroutine H_2lvl(H, omega, phi, delta)
    implicit none

    complex, dimension(2, 2), intent(out) :: H
    real, intent(in)                      :: omega, phi, delta
    complex, parameter                    :: i = (0, 1)

    H(1, 1) = 0
    H(1, 2) = omega * exp(i * phi)
    H(2, 1) = omega * exp(-i * phi)
    H(2, 2) = 2 * delta
    H = 1. / 2. * H
    
end subroutine H_2lvl

subroutine H_3lvl(H, omega_p, phi_p, delta_p, omega_s, phi_s, delta_s)
    implicit none

    complex, dimension(3, 3), intent(out) :: H
    real, intent(in)                      :: omega_p, phi_p, delta_p
    real, intent(in)                      :: omega_s, phi_s, delta_s
    complex, parameter                    :: i = (0, 1)

    H(1, 1) = 0
    H(1, 2) = omega_p * exp(i * phi_p)
    H(1, 3) = 0
    H(2, 1) = omega_p * exp(-i * phi_p)
    H(2, 2) = 2 * delta_p
    H(2, 3) = omega_s * exp(i * phi_s)
    H(3, 1) = 0
    H(3, 2) = omega_s * exp(-i * phi_s)
    H(3, 3) = 2 * (delta_p - delta_s)
    H = 1. / 2. * H
    
end subroutine H_3lvl
