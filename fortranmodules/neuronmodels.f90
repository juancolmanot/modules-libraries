module neuronmodels
    use precision
    implicit none

    contains

    function Izhikevich(u,v,a,b,c,d,I)
        use precision
        implicit none

        real(dp), intent(in)    :: u,v,a,b,c,d,I
        real(dp)                :: u1,v1,Izhikevich(2)
        
        if (v < 30._dp) then
            v1 = min(0.04_dp*v*v + 6._dp*v + 140._dp - I - u, 30._dp)
            u1 = u + a*(b*v - u)
        else if (v >= 30._dp) then
            v1 = c
            u1 = u + d
        endif
    
        Izhikevich = [u1, v1]

    end function Izhikevich


    function Rulkov(x,y,alfa,sigma,mu,I)
        use precision
        implicit none

        real(dp), intent(in)    :: x,y,alfa,sigma,mu,I
        real(dp)                :: x1,y1,Rulkov(2)

        if (x <= 0._dp) then
            x1 = alfa*(1 - x)**(-1) + y + I
        else if (x > 0._dp .and. x < alfa + y) then
            x1 = alfa + y + I
        else if (x >= alfa + y) then
            x1 = -1._dp
        endif

        y1 = y - mu*(x - sigma)

        Rulkov = [x1, y1]

    end function Rulkov

    function SupRulkov(x,y,alfa,sigma,mu,I)
        use precision
        implicit none

        real(dp), intent(in)    :: x,y,alfa,sigma,mu,I
        real(dp)                :: x1,y1,SupRulkov(2)

        if (x < -1._dp - alfa*0.5_dp) then
            x1 = -alfa*alfa*0.25_dp - alfa + y + I
        else if (x >= -1._dp - alfa*0.5_dp .and. x <= 0) then
            x1 = alfa*x + (x + 1)**2 + y + I
        else if (x > 0 .and. x < 1 + y) then
            x1 = 1._dp + y + I
        else if (x >= 1._dp + y) then
            x1 = -1._dp
        endif

        y1 = y - mu*(x - sigma)

        SupRulkov = [x1, y1]

    end function SupRulkov

    function CNV(x,y,m0,m1,a,beta,eps,d,J)
        use precision
        use functions
        implicit none

        real(dp), intent(in)    :: x,y,m0,m1,a,beta,eps,d,J
        real(dp)                :: Jmin,Jmax,x1,y1,CNV(2),F

        Jmin = a*m1/(m0 + m1); Jmax = (m0 + a*m1)/(m0 + m1)

        if (x <= Jmin) then
            F = -m0*x
        else if (x > jmin .and. x < Jmax) then
            F = m1*(x - a)
        else if (x >= Jmax) then
            F = -m0*(x - 1._dp)
        endif

        x1 = x + F - y - beta*Heaviside(x,d)
        y1 = y + eps*(x - J)

        CNV = [x1, y1]

    end function CNV

    function Chialvo(x,y,a,b,c,I)
        use precision
        implicit none
        
        real(dp), intent(in)    :: x,y,a,b,c,I
        real(dp)                :: x1,y1,Chialvo(2)

        x1 = x*x*exp(y-x) + I
        y1 = a*y - b*x + c

        Chialvo = [x1, y1]
        

    end function Chialvo

end module neuronmodels