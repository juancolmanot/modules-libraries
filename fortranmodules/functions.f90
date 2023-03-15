module functions
    use precision
    implicit none

    contains

    function Heaviside(x,d)
        use precision
        implicit none

        real(dp), intent(in)    :: x,d
        real(dp)                :: Heaviside

        if (x < d) then
            Heaviside = 0._dp
        else if (x == d) then
            Heaviside = 0.5_dp
        else if (x > d) then
            Heaviside = 1._dp
        endif 

    end function Heaviside



end module functions