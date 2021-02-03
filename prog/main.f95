program pool
implicit none
    complex*8              :: H2(2, 2), H3(3, 3)
    integer                :: i
    real*8, allocatable    :: e0(:)
    complex*8, allocatable :: phi0(:)
    logical                :: conv

    H2 = H_2lvl(1d0, 1d0, 1d0)
    call main(H2, size(H2, 1))

    H3 = H_3lvl(1d0, 1d0, 1d0, 1d0, 1d0, 1d0)
    call main(H3, size(H3, 1))

contains

    subroutine main(H, dim)
        implicit none
        integer                 :: dim
        complex, dimension(dim) :: H

        allocate(phi0(size(H, 1)))
        allocate(e0(size(H, 1)))

        call grounds_state(size(H, 1), H, 1d0, 1d-8, 1d4, phi0, e0, conv)

        do i=1, size(phi0)
            write(*, *) phi0(:)
        enddo

        print*, e0
        print*, conv

        deallocate(phi0)
        deallocate(e0)

        
    end subroutine main

    function H_2lvl(omega, phi, delta)
    implicit none

        complex*8, dimension(2, 2) :: H_2lvl
        real*8, intent(in)         :: omega, phi, delta
        complex*8, parameter       :: i = (0, 1)

        H_2lvl(1, 1) = 0d0
        H_2lvl(1, 2) = omega * exp(i * phi)
        H_2lvl(2, 1) = omega * exp(-i * phi)
        H_2lvl(2, 2) = 2d0 * delta
        H_2lvl = 1d0 / 2d0 * H_2lvl
        
    end function H_2lvl

    function H_3lvl(omega_p, phi_p, delta_p, omega_s, phi_s, delta_s)
    implicit none

        complex*8, dimension(3, 3) :: H_3lvl
        real*8, intent(in)         :: omega_p, phi_p, delta_p
        real*8, intent(in)         :: omega_s, phi_s, delta_s
        complex*8, parameter       :: i = (0d0, 1d0)

        H_3lvl(1, 1) = 0d0
        H_3lvl(1, 2) = omega_p * exp(i * phi_p)
        H_3lvl(1, 3) = 0d0
        H_3lvl(2, 1) = omega_p * exp(-i * phi_p)
        H_3lvl(2, 2) = 2d0 * delta_p
        H_3lvl(2, 3) = omega_s * exp(i * phi_s)
        H_3lvl(3, 1) = 0d0
        H_3lvl(3, 2) = omega_s * exp(-i * phi_s)
        H_3lvl(3, 3) = 2d0 * (delta_p - delta_s)
        H_3lvl = 1d0 / 2d0 * H_3lvl
        
    end function H_3lvl

end program pool
