subroutine lanczos(H, N)
    implicit none
    integer(kind=4) :: N, i
    complex(kind=8) :: H(N, N), Hb(N, N)
    complex(kind=8) :: lambda, lambda_0

    Hb = krilov(H, N)

    lambda = sum( (/ (real(Hb(i,i)), i=1, N) /) ) + 1d0
    ! p = caract_polynomiale(Hb, N, lambda)
    lambda_0 = dichotomy(Hb, N)

    do i=1, N
        print*, Hb(i, :)
    enddo

    print*, " "
    print*, lambda_0

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

    function caract_polynomiale(H, N, lambda) result(p)
        implicit none
        integer(kind=4), intent(in) :: N
        complex(kind=8), intent(in) :: H(N, N)
        real(kind=8)                :: lambda
        complex(kind=8)             :: p(N + 1)
        integer                     :: i

        p(1) = 1
        p(2) = H(1, 1) - lambda

        do i=3, N + 1
            p(i) = (H(i - 1, i - 1) - lambda) * p(i - 1) - abs(H(i - 2, i - 1))**2 * p(i - 2)
        enddo
        
    end function caract_polynomiale

    function w(x, H, N) result(count)
        implicit none
        integer(kind=4), intent(in) :: N
        real(kind=8), intent(in)    :: x
        complex(kind=8)             :: p(N + 1), H(N, N)
        integer(kind=4)             :: i, count

        p = caract_polynomiale(H, N, x)
        count = 0

        do i=1, N
            if (real(p(i + 1) / p(i)) < 0.d0) then
                count = count + 1
            endif
        enddo
        
    end function w

    function dichotomy(H, N) result(lambda_0)
        implicit none
        integer(kind=4), intent(in) :: N
        complex(kind=8), intent(in) :: H(N, N)
        real(kind=8)                :: a, b, lambda_0
        real(kind=8), parameter     :: eps = 1.d-8
        integer, parameter          :: max_iter = 1000
        integer                     :: iter

        iter = 0

        a = -1.d4
        b = 0.d0

        do while (w(b, H, N) - w(a, H, N) /= 1)
            b = (a + b) / 2.d0
            iter = iter + 1
            if (iter == max_iter) then
                print*, "The 'a' value is not negative or enought wide &
                                    &or the matrix doesn't have a negative spectrum"
                stop
            endif
        enddo

        do while (b - a > eps)
            if (w((b + a) / 2.d0, H, N) - w(a, H, N) == 1) then
                b = (b + a) / 2.d0
            else
                a = (b + a) / 2.d0
            endif
            print*, (b + a) / 2.d0
        enddo

        lambda_0 = (b + a) / 2.d0

    end function dichotomy

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