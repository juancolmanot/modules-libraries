module vector_statistics
    use precision
    use matrixalgebra
    implicit none

    contains

    function x_n_avg(x,n)
        use precision
        use matrixalgebra
        implicit none

        real(dp), intent(in)    :: x(:)
        integer(dp), intent(in) :: n
        integer(dp)             :: i,j,size_x,n_bin
        real(dp)                :: x_n_avg(n),x_n_domain(n),a
        integer(dp)             :: counter(n)
        real(dp), allocatable   :: x_sorted(:)

        size_x = size(x)

        allocate(x_sorted(size_x))
        x_sorted = 0._dp

        x_sorted = sort_minor(x)

        counter = 0_dp
        x_n_avg = 0._dp
        x_n_domain = 0._dp

        x_n_domain = linspace(x_sorted(1),x_sorted(size_x),n)

        n_bin = 1_dp

        do i = 1,size_x
            if (x_sorted(i) <= x_n_domain(n_bin+1)) then
                x_n_avg(n_bin) = x_n_avg(n_bin) + x_sorted(i)
                counter(n_bin) = counter(n_bin) + 1_dp
            else if (x_sorted(i) > x_n_domain(n_bin+1)) then
                n_bin = n_bin + 1_dp
            endif
        enddo

        do i = 1,n
            a = x_n_avg(i)
            x_n_avg(i) = x_n_avg(i)/real(counter(i),dp)
            if (isnan(x_n_avg(i)) .and. i < n) then
                x_n_avg(i) = 0._dp
            endif
        enddo
    
    end function x_n_avg


end module vector_statistics