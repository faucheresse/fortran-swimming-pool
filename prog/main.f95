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

subroutine grounds_state(dimH, H, mu, eps, kmax, phi0, e0, conv)
    implicit none
    integer           , intent(in)  :: dimH
    real, intent(in)                :: H(dimH, dimH), mu, eps, kmax
    real, intent(out)               :: phi0(dimH), e0
    logical, intent(out)            :: conv
    integer                         :: k, id(dimH, dimH), i, j
    real                            :: vect(dimH)

    conv = .false.
    id = identity(dimH)

    H = H - mu * id
    vect = prod_mat_vect(H, phi0, dimH) - dot_product(dot_product(prod_mat_vect(phi0, H, dimH), phi0), phi0)
    do while (norm(vect, dimH) > eps .and. k <= kmax)
        phi0 = prod_mat_vect(H, phi0, dimH)
        phi0 = phi0 / norm(phi0, dimH)
        k = k + 1
    enddo

    if (k > kmax) then
        conv = .true.
    endif
    H = H + mu * id
    e0 = dot_product(prod_mat_vect(phi0, H, dimH), phi0)

    contains

    function identity(dim)
        implicit none
        integer                      :: dim
        integer, dimension(dim, dim) :: identity

        do i=1, dimH
            do j=1, dimH
                if (i == j) then
                    identity(i, j) = 1
                else
                    identity(i, j) = 0
                endif
            enddo
        enddo
        
    end function identity

    real function norm(v, dimv)
        implicit none
        integer               :: dimv
        real, dimension(dimv) :: v

        do i=1, dimv
            norm = norm + v(i) * v(i)
        enddo

        norm = sqrt(norm)
        
    end function norm

    function prod_mat_vect(M, v, dim)
        implicit none
        integer :: dim
        real    :: M(dim, dim), v(dim), prod_mat_vect(dim)

        do i=1, dim
            prod_mat_vect(i) = dot_product(M(i, :), v)
        enddo
        
    end function prod_mat_vect
    
end subroutine grounds_state

