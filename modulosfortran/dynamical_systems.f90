module dynamical_systems
    use precision
    implicit none

    contains

    function map_A(x,eps,a,d)
        implicit none

        real(dp), intent(in)    :: x,eps,a,d
        real(dp)                :: map_A,x3,x6

        x3 = x*x*x; x6 = x3*x3

        map_A = -(1._dp + eps)*x - a*x3 + d*x6*sin(x)

    end function   map_A




end module dynamical_systems