program pool
implicit none
    complex(kind=8) :: H2(2, 2)
    ! complex(kind=8) :: H3(3, 3)
    

    H2 = H_2lvl(2.d0, 0.d0, 1.d0)
    call main(H2, size(H2, 1))

    ! H3 = H_3lvl(0.d0, 8.d-1, 5d-1, 0.d-8, 8.d-1, -5d-1)
    ! call main(H3, size(H3, 1))

contains

    subroutine main(H, dim)
        implicit none
        integer                      :: dim, i
        complex(kind=8)              :: H(dim, dim)
        real(kind=8)                 :: shift
        real(kind=8)                 :: eig_0, eig_1
        complex(kind=8), allocatable :: phi0(:), phi1(:)
        logical                      :: conv_0, conv_1

        allocate(phi0(size(H, 1)))
        allocate(phi1(size(H, 1)))

        ! do i=1, dim
        !     print*, H(i, :)
        ! enddo

        shift = sum( (/ (real(H(i,i)), i=1, size(H, 1)) /) ) + 1d0
        phi0 = phi0 * 0d0
        phi0(1) = (1.d0, 0d0)
        phi1 = phi1 * 0d0
        phi1(2) = (1.d0, 0d0)

        call state(size(H, 1), H, shift, 1d-8, 1d3, phi0, eig_0, conv_0, phi1, eig_1, conv_1)

        deallocate(phi0)
        deallocate(phi1)

        
    end subroutine main

    function H_2lvl(omega, phi, delta)
    implicit none

        complex(kind=8), dimension(2, 2) :: H_2lvl
        real(kind=8), intent(in)         :: omega, phi, delta
        complex(kind=8), parameter       :: i = (0, 1)

        H_2lvl(1, 1) = 0d0
        H_2lvl(1, 2) = omega * exp(i * phi)
        H_2lvl(2, 1) = omega * exp(-i * phi)
        H_2lvl(2, 2) = 2d0 * delta
        H_2lvl = 1d0 / 2d0 * H_2lvl
        
    end function H_2lvl

    function H_3lvl(omega_p, phi_p, delta_p, omega_s, phi_s, delta_s)
    implicit none

        complex(kind=8), dimension(3, 3) :: H_3lvl
        real(kind=8), intent(in)         :: omega_p, phi_p, delta_p
        real(kind=8), intent(in)         :: omega_s, phi_s, delta_s
        complex(kind=8), parameter       :: i = (0, 8)

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

    subroutine molecular_potential(N, L, V0, x0, a, V)
        implicit none
        integer, intent(in) :: N
        real, intent(in)    :: L
        real, intent(out)   :: V(N, N)
        integer             :: i
        real                :: x(N)

        V = V * 0.d0

        x(1) = 0
        do i=2, N
            x(i) = x(i - 1) + L / N
        enddo

        ! do i=1, N
        !     V(i, i) = 
        ! enddo
        
    end subroutine molecular_potential

    ! subroutine molecular_H(L)
    !     implicit none
        
    ! end subroutine molecular_H

end program pool
