subroutine lanczos(H, N)
    implicit none
    integer(kind=4) :: N, i
    complex(kind=8) :: H(N, N), Hb(N, N)

    Hb = krilov(H, N)

    do i=1, N
        print*, Hb(i, :)
    enddo

contains

    function krilov(H, N) result(Hb)
        implicit none
        integer(kind=4), intent(in) :: N
        complex(kind=8), intent(in) :: H(N, N)
        complex(kind=8)             :: Hb(N, N)
        complex(kind=8)             :: psi(N, N + 1), alpha, beta, s(N)
        integer                     :: i, j

        Hb = Hb * 0.d0
        s = s * 0.d0

        psi = psi * 0.d0
        psi(1, 1) = 1.d0
        

        alpha = dot_product(psi(:, 1), matmul(H, psi(:, 1)))
        Hb(1, 1) = alpha
        psi(:, 2) = matmul(H, psi(:, 1)) - alpha * psi(:, 1)
        psi(:, 2) = psi(:, 2) / norm(psi(:, 2), N)

        do i=1, N - 1
            alpha = dot_product(psi(:, i + 1), matmul(H, psi(:, i + 1)))
            Hb(i + 1, i + 1) = alpha
            do j=1, i + 2
                s = s + psi(:, j)
            enddo
            psi(:, i + 2) = matmul(H, psi(:, i + 1)) - alpha * s
            psi(:, i + 2) = psi(:, i + 2) / norm(psi(:, i + 2), N)

            beta  = dot_product(psi(:, i + 2), matmul(H, psi(:, i + 1)))
            Hb(i, i + 1) = conjg(beta)
            Hb(i + 1, i) = beta
        enddo
        
    end function krilov

    real(kind=8) function norm(v, dimv)
        implicit none
        integer(kind=4) :: dimv
        complex(kind=8) :: v(dimv)

        norm = sqrt(real(dot_product(v, v)))

        if (norm == 0.) then
            norm = 1.
        endif
        
    end function norm

end subroutine lanczos