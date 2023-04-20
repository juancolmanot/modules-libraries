module maps
    use precision
    implicit none
    
    contains

    ! Type-I intermittency maps

    ! Elaskar - Del Río map from book (page 11)
    function map_I_a(x,eps,a)
        use precision
        implicit none

        real(dp), intent(in)    :: x,eps,a
        real(dp)                :: map_I_a

        map_I_a = eps + x + a*x*x

    end function map_I_a

    ! Type-II intermittency maps

    ! Elaskar - Del Río map from book (page 19)
    function map_II_a(x,eps,a,b,q)
        use precision
        implicit none

        real(dp), intent(in)    :: x(2),eps,a,b,q
        real(dp)                :: map_II_a(2),x3

        x3 = x(1)*x(1)*x(1) 

        map_II_a(1) = (1._dp + eps)*x(1) + a*x3
        map_II_a(1) = x(2) + b + q*x(1)*x(1)

    end function map_II_a

    ! Type-III intermittency maps

    ! Elaskar - Del Río map from book (page 15)
    function map_III_a(x,eps,a,d)
        use precision
        implicit none

        real(dp), intent(in)    :: x,eps,a,d
        real(dp)                :: map_III_a,x3,x6

        x3 = x*x*x; x6 = x3*x3

        map_III_a = -(1._dp + eps)*x - a*x3 + d*x6*sin(x)

    end function map_III_a

    
end module maps