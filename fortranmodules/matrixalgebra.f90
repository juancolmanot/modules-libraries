module matrixalgebra
    use precision
    implicit none

    contains

    function linspace(a,b,n)
        use precision
        implicit none

        integer(dp), intent(in) :: n
        real(dp), intent(in)    :: a,b
        integer(sp)             :: i
        real(dp)                :: linspace(n)

        linspace = 0._dp

        do i = 1,n
            linspace(i) = a + (b - a)*real(i-1,dp)/real(n-1,dp)
        enddo
    
    end function linspace

    function sort(x)
        use precision
        implicit none

        real(dp), intent(in)    :: x(:)
        integer(dp)             :: i,n,min_idx
        real(dp), allocatable   :: sort(:),x_1(:)
        real(dp)                :: x_aux

        n = size(x)

        allocate(sort(n),x_1(n))

        x_1 = x

        do i = 1,n
            sort(i) = minval(x_1(i:n))
            min_idx = minloc(x_1(i:n), 1)
            x_aux = x_1(i)
            x_1(i) = sort(i)
            x_1(i - 1 + min_idx) = x_aux
        enddo

    end function sort

end module matrixalgebra
