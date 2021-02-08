subroutine state(dimH, H, mu, eps, kmax, phi0, eig_0, conv_0, phi1, eig_1, conv_1)
    implicit none
    integer*4, intent(in)    :: dimH
    complex*8                :: H(dimH, dimH)
    real*8, intent(in)       :: mu, eps, kmax
    complex*8, intent(inout) :: phi0(dimH), phi1(dimH)
    real*8, intent(out)      :: eig_0, eig_1
    logical, intent(out)     :: conv_0, conv_1
    integer*8                :: k, id(dimH, dimH)

    call ground_state
    call first_state

    print"('case: (', I1, ', ', I1, ')')", size(H, 1), size(H, 2)
    print*, "phi0:", phi0
    print*, "phi1:", phi1
    print*, "eig_0 = ", eig_0
    print*, "eig_1 = ", eig_1
    print*, "Has converged: ground state: ", conv_0, " first_state: ", conv_1
    print*, " "

    contains

    subroutine ground_state
        implicit none

        conv_0 = .true.
        id = identity(dimH)
        k = 0

        phi0 = phi0 / norm(phi0, dimH)
        H = H - mu * id

        do while (norm(matmul(H, phi0) - dot_product(phi0, matmul(H, phi0))&
                                            * phi0, dimH) > eps .and. k < kmax)
            phi0 = matmul(H, phi0)
            phi0 = phi0 / norm(phi0, dimH)
            k = k + 1

        enddo

        if (k >= kmax) then
            conv_0 = .false.
        endif

        H = H + mu * id

        eig_0 = real(dot_product(phi0, matmul(H, phi0)))
        
    end subroutine ground_state

    subroutine first_state
        implicit none

        conv_1 = .true.
        id = identity(dimH)
        k = 0

        phi1 = phi1 - dot_product(phi0, phi1) * phi0
        phi1 = phi1 / norm(phi1, dimH)

        H = H - mu * id

        do while (norm(matmul(H, phi1) - dot_product(phi1, matmul(H, phi1))&
                                            * phi1, dimH) > eps .and. k < kmax)
            phi1 = matmul(H, phi1)
            phi1 = phi1 - dot_product(phi0, phi1) * phi0
            phi1 = phi1 / norm(phi1, dimH)
            k = k + 1

        enddo

        if (k >= kmax) then
            conv_1 = .false.
        endif

        H = H + mu * id

        eig_1 = real(dot_product(phi1, matmul(H, phi1)))
        
    end subroutine first_state

    function identity(dim)
        implicit none
        integer*4                      :: dim, i, j
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
        integer*4  :: dimv
        complex*8  :: v(dimv)

        norm = sqrt(real(dot_product(v, v)))
        
    end function norm
    
end subroutine state
