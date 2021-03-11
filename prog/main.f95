program pool
implicit none
    ! complex(kind=8)              :: H2(2, 2)
    integer :: i
    complex(kind=8)              :: H3(3, 3)
    ! complex(kind=8), allocatable :: mol_H(:, :)
    ! integer, parameter           :: N = 3 !has to be positive
    ! complex(kind=8) :: test(3, 3)
    
    ! test(1, 1) = 51.5
    ! test(1, 2) = 105.5
    ! test(1, 3) = 54.5
    ! test(2, 1) = -15. 
    ! test(2, 2) = -32.
    ! test(2, 3) = -15.
    ! test(3, 1) = -28.5
    ! test(3, 2) = -55.5
    ! test(3, 3) = -31.5

    ! print*, "Power iterative method :"
    ! call main(test, size(test, 1))
    ! print*, "Lanczos Algorithm :"
    ! call lanczos(test, size(test, 1))
    

    ! print*, ""
    ! H2 = H_2lvl(i/10000.d0, 0.d0, 1.d0)
    ! print*, "Power iterative method :"
    ! call main(H2, size(H2, 1), i/10000.d0)
    ! print*, "Lanczos Algorithm :"
    ! call lanczos(H2, size(H2, 1))

    ! print*, ""

    ! H3 = H_3lvl(0.d0, 8.d-1, 5d-1, 0.d0, 8.d-1, -5d-1)
    ! print*, "Power iterative method :"
    ! call main(H3, size(H3, 1))
    ! print*, "Lanczos Algorithm :"
    ! call lanczos(H3, size(H3, 1))

    ! print*, ""

    ! allocate(mol_H(N, N))
    ! mol_H = molecular_H(N, 10.d0, 2.01588d0)
    ! mol_H = mol_H + molecular_potential(N, 10.d0, 5.d0, 2.d0)
    ! ! mol_H = mol_H + box_potential(N, 10.d0)
    ! ! mol_H = mol_H + ho_potential(N, 10.d0, -5.d-1)
    ! print*, "Power iterative method :"
    ! call main(mol_H, N)
    ! ! print*, "Lanczos Algorithm :"
    ! ! call lanczos(mol_H, N)
    ! deallocate(mol_H)

contains

    subroutine main(H, dim, om)
        implicit none
        integer                      :: dim, i
        complex(kind=8)              :: H(dim, dim)
        real(kind=8)                 :: shift
        real(kind=8)                 :: eig_0, eig_1, om
        complex(kind=8), allocatable :: phi0(:), phi1(:)
        logical                      :: conv_0, conv_1

        allocate(phi0(size(H, 1)))
        allocate(phi1(size(H, 1)))

        ! do i=1, dim
        !     print*, H(i, :)
        ! enddo
        ! print*, ""

        shift = sum( (/ (real(H(i,i)), i=1, size(H, 1)) /) ) + 1d0
        phi0 = phi0 * 0d0
        phi0(1) = (1.d0, 0d0)
        phi1 = phi1 * 0d0
        phi1(2) = (1.d0, 0d0)

        call state(size(H, 1), H, shift, 1d-8, 1d3, phi0, eig_0, conv_0, phi1, eig_1, conv_1)

        open(unit=10, file="data_H3.dat",position="append")

        write(10, *) om, real(phi0)
        close(10)



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

    function molecular_potential(N, L, V0, a)  result(V)
        implicit none
        integer, intent(in)          :: N
        real(kind=8), intent(in)     :: L, V0, a
        complex(kind=8)              :: V(N, N)
        integer                      :: i
        real(kind=8)                 :: x(N)

        V = V * 0.d0

        do i=1, N
            x(i) = (i - 1) * L / (N - 1)
        enddo

        do i=1, N
            V(i, i) = V0 * (exp(-2 * a * x(i)) - 2 * exp(-a * x(i)))
        enddo
        
    end function molecular_potential

    function box_potential(N, L) result(V)
        implicit none
        integer, intent(in)          :: N
        real(kind=8), intent(in)     :: L
        complex(kind=8)              :: V(N, N)
        integer                      :: i
        real(kind=8)                 :: x(N)

        V = V * 0.d0

        do i=1, N
            x(i) = (i - 1) * L / (N - 1)
        enddo

        do i=1, N
            if (x(i) == 0. .or. x(i) == L) then
                v(i, i) = 0.
            else
                v(i, i) = 10.d10
            endif
        enddo
        
    end function box_potential

    function ho_potential(N, L, k) result(V)
        implicit none
        integer, intent(in)          :: N
        real(kind=8), intent(in)     :: L, k
        complex(kind=8)              :: V(N, N)
        integer                      :: i
        real(kind=8)                 :: x(N)

        V = V * 0.d0

        do i=1, N
            x(i) = (i - 1) * L / (N - 1)
        enddo

        do i=1, N
            v(i, i) = 1. / 2. * k * (x(i) - L / 2.)**2
        enddo
        
    end function ho_potential

    function molecular_H(N, L, m) result(H)
        implicit none
        integer(kind=4), intent(in)      :: N
        real(kind=8), intent(in)         :: L, m
        complex(kind=8), dimension(N, N) :: H
        real(kind=8), parameter          :: pi = 4 * atan(1.d0)
        integer                          :: i, j

        do i=1, N
            do j=1, N
                if (i .eq. j) then
                    H(i, j) = (pi**2. * (N + 1.)**2. + 2.) / 3. / L**2.
                else
                    H(i, j) = ((-1.)**(j - i) * 2. * pi**2.) / (L**2. * sin((j - i) * pi / (N + 1.))**2.)
                endif
            enddo
        enddo

        H = H / 2. / m
        
    end function molecular_H

end program pool
