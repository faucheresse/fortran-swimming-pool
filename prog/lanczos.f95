subroutine lanczos(H, N)
    implicit none
    integer(kind=4) :: N, i
    complex(kind=8) :: H(N, N), Hb(N, N)
    complex(kind=8) :: lambda_0

    call krilov
    call dichotomy

    do i=1, N
        print*, Hb(i, :)
    enddo

    print*, " "
    print*, lambda_0

contains

    subroutine krilov
        implicit none
        complex(kind=8)             :: psi(N, N + 1), s(N)
        integer                     :: i, k

        Hb = Hb * 0.d0

        psi = psi * 0.d0
        psi(1, 1) = 1.d0

        Hb(1, 1) = dot_product(psi(:, 1), matmul(H, psi(:, 1)))

        psi(:, 2) = matmul(H, psi(:, 1)) - dot_product(psi(:, 1), matmul(H, psi(:, 1))) * psi(:, 1)
        psi(:, 2) = psi(:, 2) / norm(psi(:, 2), N)

        do i=1, N - 1
            s = s * 0.d0
            do k=1, i
                s = s + dot_product(psi(:, k), matmul(H, psi(:, i + 1))) * psi(:, k)
            enddo

            psi(:, i + 2) = matmul(H, psi(:, i + 1)) - s
            psi(:, i + 2) = psi(:, i + 2) / norm(psi(:, i + 2), N)

            Hb(i + 1, i + 1) = dot_product(psi(:, i + 1), matmul(H, psi(:, i + 1)))
            Hb(i + 1, i) = dot_product(psi(:, i + 2), matmul(H, psi(:, i + 1)))
            Hb(i, i + 1) = conjg(Hb(i + 1, i))
        enddo

    end subroutine krilov

    subroutine dichotomy
        implicit none
        real(kind=8)                 :: a, b
        real(kind=8), parameter      :: eps = 1.d-8
        integer, parameter           :: max_iter = 1000
        integer                      :: iter

        iter = 0

        a = -1.d8
        b = 0.d0

        do while (w(b) - w(a) /= 1)
            b = (a + b) / 2.d0
            iter = iter + 1
            if (iter == max_iter) then
                print*, "Max iteration reached"
                exit
            endif
        enddo

        do while (b - a > eps)
            if (w((b + a) / 2.d0) - w(a) == 1) then
                b = (b + a) / 2.d0
            else
                a = (b + a) / 2.d0
            endif
            ! print*, (b + a) / 2.d0
        enddo

        lambda_0 = (b + a) / 2.d0

    end subroutine dichotomy

    function caract_polynomial(lambda) result(p)
        implicit none
        real(kind=8), intent(in) :: lambda
        complex(kind=8)          :: p(N + 1)
        integer                  :: i

        p(1) = 1
        p(2) = H(1, 1) - lambda

        do i=3, N + 1
            p(i) = (Hb(i - 1, i - 1) - lambda) * p(i - 1) - abs(Hb(i - 2, i - 1))**2 * p(i - 2)
        enddo
        
    end function caract_polynomial

    function w(x) result(count)
        implicit none
        real(kind=8), intent(in)    :: x
        complex(kind=8)             :: p(N + 1)
        integer(kind=4)             :: i, count

        p = caract_polynomial(x)
        count = 0

        do i=1, N
            if (real(p(i + 1) / p(i)) < 0.d0) then
                count = count + 1
            endif
        enddo
        
    end function w

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