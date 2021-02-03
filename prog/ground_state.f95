subroutine grounds_state(dimH, H, mu, eps, kmax, phi0, e0, conv)
    implicit none
    integer, intent(in)  :: dimH
    complex*8            :: H(dimH, dimH), phi0(dimH)
    real*8, intent(in)   :: mu, eps, kmax
    real*8, intent(out)  :: e0(dimH)
    logical, intent(out) :: conv
    integer              :: k, id(dimH, dimH), i, j
    complex*8            :: vect(dimH)

    conv = .false.
    id = identity(dimH)

    H = H - mu * id
    vect = mat_dot_product(H, phi0, dimH) &
            - mat_dot_product(trans_dot_product(phi0,&
                mat_dot_product(H, phi0, dimH), dimH), phi0, dimH)
    do while (norm(vect, dimH) > eps .and. k <= kmax)
        phi0 = mat_dot_product(H, phi0, dimH)
        phi0 = phi0 / norm(phi0, dimH)
        k = k + 1

    enddo

    if (k > kmax) then
        conv = .true.
    endif
    H = H + mu * id
    e0 = dot_product(phi0, mat_dot_product(H, phi0, dimH))

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

    real*8 function norm(v, dimv)
        implicit none
        integer :: dimv
        complex*8  :: v(dimv)

        do i=1, dimv
            norm = norm + conjg(v(i)) * v(i)
        enddo

        norm = sqrt(norm)
        
    end function norm

    function mat_dot_product(M, v, dim)
        implicit none
        integer :: dim, i
        complex*8 :: M(dim, dim), mat_dot_product(dim), v(dim)

        do i=1, dim
            mat_dot_product(i) = dot_product(M(i, :), v)
        enddo
        
    end function mat_dot_product

    function trans_dot_product(v1, v2, dim)
        implicit none
        integer :: dim, i
        complex*8 :: v1(dim), v2(dim),trans_dot_product(dim)

        do i=1, dim
            trans_dot_product(i) = v1(i) * v2(i)
        enddo
        
    end function trans_dot_product
    
end subroutine grounds_state
