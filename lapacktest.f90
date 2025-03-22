program lapack_test
    implicit none
    integer, parameter :: n = 3
    real(8) :: a(n, n), b(n), x(n)
    integer :: ipiv(n), info, i

    ! Example matrix and vector
    a = reshape((/3.0, 1.0, 2.0, 1.0, 4.0, 1.0, 2.0, 1.0, 5.0/), (/n, n/))
    b = (/1.0, 2.0, 3.0/)

    ! Solve the linear system Ax = b using DGESV
    call dgesv(n, 1, a, n, ipiv, b, n, info)

    if (info == 0) then
        print *, "Solution:"
        do i = 1, n
            print *, b(i)
        end do
    else
        print *, "Error solving linear system, info =", info
    end if

end program lapack_test