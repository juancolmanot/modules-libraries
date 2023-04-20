module intermittency_type_III_DR_E_D_2011
    use precision
    implicit none

    contains

    !=======================================================
    ! Type III intermittency functions
    !-------------------------------------------------------
    ! M function theoric only positives
    function Mx(x, m, xi, c)
        use precision
        implicit none
        
        real(dp), intent(in) :: x, m, xi, c
        real(dp)             :: Mx

        if (x < xi .or. x > c) then
            Mx = 0._dp
        else
            Mx = m * (x - xi) + xi
        endif

    end function Mx
    
    !-------------------------------------------------------
    ! RPD one side
    function RPD(x, xi, c, m)
        use precision
        implicit none

        real(dp), intent(in)    :: x, xi, c, m
        real(dp)                :: RPD, b, alfa

        alfa = (1._dp - 2._dp * m) / (m - 1._dp)

        b = 0.5_dp * (alfa + 1._dp) / ((c - xi)**(int(alfa, dp) + 1_dp))

        if (x < xi .or. x > c) then
            RPD = 0._dp
        else
            RPD = b * (x - xi)**(int(alfa, dp))
        endif

    end function RPD

    !-------------------------------------------------------
    ! RPD absolute
    function RPD_abs(x, xi, c, alfa)
        use precision
        implicit none

        real(dp), intent(in)    :: x, xi, c, alfa
        real(dp)                :: RPD_abs, b

        b = 0.5_dp * (alfa + 1._dp) / ((c + abs(xi))**(int(alfa, dp) + 1_dp))

        if (abs(x) < abs(xi)) then
            RPD_abs = b * ((abs(xi) + x)**int(alfa, dp) + (abs(xi) - x)**(int(alfa, dp)))
        else if (abs(xi) < x .and. x <= c) then
            RPD_abs = b * (abs(xi) + x)**int(alfa, dp)
        else if (-c < x .and. x <= -abs(xi)) then
            RPD_abs = b * (abs(xi) - x)**int(alfa, dp)
        endif

    end function RPD_abs

    !-------------------------------------------------------
    ! laminar length
    function laminarlength(x, c, a, eps)
        use precision
        implicit none

        real(dp), intent(in)    :: x, c, a, eps
        real(dp)                :: laminarlength

        laminarlength = 0.5_dp * ((2._dp * log(c / abs(x)) - log((eps + a * c**2) / (eps + a * x**2)))) / eps

    end function laminarlength
    
    !-------------------------------------------------------
    ! Inverse X(l, c)
    function inverse_X(l, c, a, eps)
        use precision
        implicit none

        real(dp), intent(in)    :: l, c, a, eps
        real(dp)                :: inverse_X

        inverse_X = sqrt(eps / ((a + (eps / c**2)) * exp(2._dp * eps * l) - a))

    end function inverse_X    

    !-------------------------------------------------------
    ! laminar length probability
    function laminar_prob(x, l, alfa, c, a, eps)
        use precision
        implicit none

        real(dp), intent(in)    :: l, alfa, c, a, eps
        real(dp)                :: laminar_prob, x

        x = inverse_X(l, c, a, eps)

        laminar_prob = 2._dp * b * ((x - xi)**int(alfa, dp)) * (a * x**3 + eps * x)

    end function laminar_prob
    !=======================================================




end module intermittency_type_III_DR_E_D_2011