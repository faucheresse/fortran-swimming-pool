subroutine lanczos(H, dim, psi)
    implicit none
    integer(kind=4) :: dim
    complex(kind=8) :: H(N, N), psi(N)

    

contains

    subroutine krilov
        implicit none
        complex(kind=8) :: psi_m(N),psi_e(N), psi_p(N)
        integer         :: i

        psi_m = psi
        psi_m = psi_m / norm(psi_m, N)
        psi_e = matmul(H, psi_m) - dot_product(psi_m, mtmul(H, psi_m)) * psi_m
        psi_e = psi_e / norm(psi_e, N)
        psi_p = matmul(H, psi_e) - dot_product(psi_e, mtmul(H, psi_e)) * psi_e
        psi_p = psi_p / norm(psi_p, N)

        do i=1, N
            H = dot_product(psi_m, mtmul(H, psi_e)) * psi_m &
                            + dot_product(psi_e, mtmul(H, psi_e)) &
                            * psi_e + dot_product(psi_p, mtmul(H, psi_e)) * psi_p


            psi_m = psi_e
            psi_m = psi_m / norm(psi_m, N)
            psi_e = psi_p
            psi_e = psi_e / norm(psi_e, N)
            psi_p = matmul(H, psi_e) - dot_product(psi_e, mtmul(H, psi_e)) * psi_e
            psi_p = psi_p / norm(psi_p, N)
        enddo
        
    end subroutine krilov

    real(kind=8) function norm(v, dimv)
        implicit none
        integer(kind=4) :: dimv
        complex(kind=8) :: v(dimv)

        norm = sqrt(real(dot_product(v, v)))
        
    end function norm

end subroutine lanczos