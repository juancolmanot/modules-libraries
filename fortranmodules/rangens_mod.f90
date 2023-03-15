module rangens
    use precision, only: dp=>dp,sp=>sp
    implicit none

    contains

    function ran0(idum)
        implicit none
        real(kind(1.d0))       :: ran0
        integer, intent(inout) :: idum
        integer                :: ia,im,iq,ir,k
        real(kind(1.d0))       :: am
        parameter(ia=16807,im=2147483647.d0,am=1.d0/2147483647.d0, &
        &   iq=127773,ir=2836)

        k = idum/iq
        idum = ia*(idum - k*iq) - ir*k
        if (idum <= 0) idum = idum + im
        ran0 = am*real(idum,kind(1.d0))
    end function ran0

    function gaussdev(idum)
        use precision
        use ran2mod, only:ran2
        implicit none

        real(dp), parameter     :: pi=4._dp*atan(1._dp)
        integer(sp), intent(in) :: idum
        real(dp)                :: u1,u2
        integer(dp)             :: iset=0
        real(dp)                :: fs,ang,gaussdev,gset
        save iset,gset

        if (iset == 0) then
            u1 = ran2(idum)
            u2 = ran2(idum)
            fs = sqrt(-2._dp*log(u1))
            ang = 2._dp*pi*u2
            gset = fs*cos(ang)
            gaussdev = fs*sin(ang)
            iset = 1
        else
            gaussdev = gset
            iset = 0
        endif
    end function gaussdev

end module rangens