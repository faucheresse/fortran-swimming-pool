program pool
implicit none
    complex*8 :: H2(2, 2)
    ! complex*8 :: H3(3, 3)
    

    H2 = H_2lvl(0.d0, 1.d0, 1d0)
    call main(H2, size(H2, 1))

    ! H3 = H_3lvl(3.5_8, 2d0 / 3d0, 5d-1, 3.5_8, 7d0 / 3d0, -5d-1)
    ! call main(H3, size(H3, 1))

contains

    subroutine main(H, dim)
        implicit none
        integer                :: dim, i
        complex*8              :: H(dim, dim)
        real*8                 :: shift
        real*8                 :: eig_0, eig_1
        complex*8, allocatable :: phi0(:), phi1(:)
        logical                :: conv_0, conv_1

        allocate(phi0(size(H, 1)))
        allocate(phi1(size(H, 1)))

        shift = sum( (/ (real(H(i,i)), i=1, size(H, 1)) /) ) + 1d0
        phi0 = phi0 * 0d0 + (2.d0, 0d0)
        phi1 = phi1 * 0d0 + (2.d0, 0d0)

        call state(size(H, 1), H, shift, 1d-8, 1d3, phi0, eig_0, conv_0, phi1, eig_1, conv_1)

        deallocate(phi0)
        deallocate(phi1)

        
    end subroutine main

    function H_2lvl(omega, phi, delta)
    implicit none

        complex*8, dimension(2, 2) :: H_2lvl
        real*8, intent(in)         :: omega, phi, delta
        complex*8, parameter       :: i = (0_8, 1_8)

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
        complex*8, parameter       :: i = (0_8, 1_8)

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
