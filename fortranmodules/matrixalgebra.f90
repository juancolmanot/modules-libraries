module matrixalgebra
    use precision
    implicit none

    contains

    function linspace(a,b,n)
        use precision
        implicit none

        integer(dp), intent(in) :: n
        real(dp), intent(in)    :: a,b
        integer(dp)             :: i
        real(dp)                :: linspace(n)

        linspace = 0._dp

        do i = 1,n
            linspace(i) = a + (b - a)*real(i-1,dp)/real(n-1,dp)
        enddo
    
    end function linspace

    function sort_minor(x)
        use precision
        implicit none

        real(dp), intent(in)    :: x(:)
        integer(dp)             :: i,n,min_idx
        real(dp), allocatable   :: sort_minor(:),x_1(:)
        real(dp)                :: x_aux

        n = size(x)

        allocate(sort_minor(n),x_1(n))

        x_1 = x

        do i = 1,n
            sort_minor(i) = minval(x_1(i:n))
            min_idx = minloc(x_1(i:n), 1)
            x_aux = x_1(i)
            x_1(i) = sort_minor(i)
            x_1(i - 1 + min_idx) = x_aux
        enddo

    end function sort_minor

    function sort_major(x)
        use precision
        implicit none

        real(dp), intent(in)    :: x(:)
        integer(dp)             :: i,n,min_idx
        real(dp), allocatable   :: sort_major(:),x_1(:)
        real(dp)                :: x_aux

        n = size(x)

        allocate(sort_major(n),x_1(n))

        x_1 = x

        do i = 1,n
            sort_major(i) = maxval(x_1(i:n))
            min_idx = maxloc(x_1(i:n), 1)
            x_aux = x_1(i)
            x_1(i) = sort_major(i)
            x_1(i - 1 + min_idx) = x_aux
        enddo

    end function sort_major

end module matrixalgebra
