program main
implicit none
    complex :: H2(2, 2), H3(3, 3)
    integer :: i

    H2 = H_2lvl(1., 1., 1.)
    H3 = H_3lvl(1., 1., 1., 1., 1., 1.)

    do i=1, size(H3, 1)
        write(*, *) H3(i, :)
    enddo

contains

    function H_2lvl(omega, phi, delta)
    implicit none

        complex, dimension(2, 2) :: H_2lvl
        real, intent(in)         :: omega, phi, delta
        complex, parameter       :: i = (0, 1)

        H_2lvl(1, 1) = 0
        H_2lvl(1, 2) = omega * exp(i * phi)
        H_2lvl(2, 1) = omega * exp(-i * phi)
        H_2lvl(2, 2) = 2 * delta
        H_2lvl = 1. / 2. * H_2lvl
        
    end function H_2lvl

    function H_3lvl(omega_p, phi_p, delta_p, omega_s, phi_s, delta_s)
    implicit none

        complex, dimension(3, 3) :: H_3lvl
        real, intent(in)         :: omega_p, phi_p, delta_p
        real, intent(in)         :: omega_s, phi_s, delta_s
        complex, parameter       :: i = (0, 1)

        H_3lvl(1, 1) = 0
        H_3lvl(1, 2) = omega_p * exp(i * phi_p)
        H_3lvl(1, 3) = 0
        H_3lvl(2, 1) = omega_p * exp(-i * phi_p)
        H_3lvl(2, 2) = 2 * delta_p
        H_3lvl(2, 3) = omega_s * exp(i * phi_s)
        H_3lvl(3, 1) = 0
        H_3lvl(3, 2) = omega_s * exp(-i * phi_s)
        H_3lvl(3, 3) = 2 * (delta_p - delta_s)
        H_3lvl = 1. / 2. * H_3lvl
        
    end function H_3lvl

end program main

