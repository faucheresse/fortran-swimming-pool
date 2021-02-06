subroutine grounds_state(dimH, H, mu, eps, kmax, phi0, eig_0, conv)
    implicit none
    integer*8, intent(in) :: dimH
    complex*8             :: H(dimH, dimH), phi0(dimH)
    real*8, intent(in)    :: mu, eps, kmax
    real*8, intent(out)   :: eig_0
    logical, intent(out)  :: conv
    integer*8             :: k, id(dimH, dimH), i, j
    complex*8             :: vect(dimH)

    conv = .false.
    id = identity(dimH)
    k = 0

    H = H - mu * id
    vect = matmul(H, phi0) - dot_product(phi0, matmul(H, phi0)) * phi0


    do while (norm(vect, dimH) > eps .and. k < kmax)
        phi0 = matmul(H, phi0)
        phi0 = phi0 / norm(phi0, dimH)
        k = k + 1
        print*, k

    enddo

    if (k >= kmax) then
        conv = .true.
    endif

    H = H + mu * id

    eig_0 = real(dot_product(phi0, matmul(H, phi0)))

    contains

    function identity(dim)
        implicit none
        integer*8                      :: dim
        integer*8, dimension(dim, dim) :: identity

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
        integer*8  :: dimv
        complex*8  :: v(dimv)

        norm = sqrt(dot_product(v, conjg(v)))
        
    end function norm
    
end subroutine grounds_state
