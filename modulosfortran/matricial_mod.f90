module matricial
    use precision, only: dp=>dp,sp=>sp
    implicit none
    
contains

subroutine tridiagonal(c,d,a,b,x)
    implicit none
    real(dp), allocatable, intent(in)   :: c(:),d(:),a(:),b(:)
    real(dp), allocatable, intent(out)  :: x(:)
    real(dp), allocatable               :: h(:),p(:)
    integer                             :: i,n

    n = size(d)

    allocate(x(n),h(n-1),p(n))

    h(1) = c(1)/d(1)
    p(1) = b(1)/d(1)
    
    do i = 2,n-1
        h(i) = c(i)/(d(i)-a(i)*h(i-1))
        p(i) = (b(i)-a(i)*p(i-1))/(d(i)-a(i)*h(i-1))
    enddo

    p(n) = (b(n)-a(n)*p(n-1))/(d(n)-a(n)*h(n-1))

    x(n) = p(n)
    
    do i = n-1,1,-1
        x(i) = p(i) - h(i)*x(i+1)  
    enddo

    deallocate(h,p)

end subroutine tridiagonal
   
end module matricial



