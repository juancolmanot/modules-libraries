module mechanics
use precision, only: pr => dp
implicit none

contains

subroutine Lagrangian(a,b,gama,tita1,tita2,omega1,omega2,w1,w2,dw1,dw2)
    use precision, only: pr => dp
    implicit none
    real(pr)                :: t1t2,o12,o22,s2t1,s2t2,s2t1t2
    real(pr), intent(in)    :: a,b,gama
    real(pr), intent(in)    :: tita1,tita2,omega1,omega2
    real(pr), intent(out)   :: w1,w2,dw1,dw2

    t1t2 = tita1-tita2
    o12 = omega1*omega1
    o22 = omega2*omega2
    s2t1 = sin(tita1)*sin(tita1)
    s2t2 = sin(tita2)*sin(tita2)
    s2t1t2 = sin(t1t2)*sin(t1t2)
    w1 = omega1
    w2 = omega2
    dw1 = (-(a+1._pr)*gama*sin(tita1)-a*b*o22*sin(t1t2)-a*cos(t1t2)*(o12*sin(t1t2)-gama*sin(tita2)))/(1._pr+a*s2t1t2)
    dw2 = ((1._pr+a)*(o12*sin(t1t2)-gama*sin(tita2))+cos(t1t2)*((1._pr+a)*gama*sin(tita1)+a*b*o22*sin(t1t2)))&
    /(b*(1._pr+a*(s2t1t2)))
    


end subroutine Lagrangian


subroutine energy_dble_pendulum(tita1,tita2,omega1,omega2,a,b,gama,T,U)
    implicit none
    real(pr), intent(in)   :: tita1,tita2,omega1,omega2,a,b,gama
    real(pr), intent(out)  :: T,U

    real(pr)               :: o12,o22,b2,t12
    
    o12 = omega1*omega1
    o22 = omega2*omega2
    t12 = tita1-tita2
    b2 = b*b

    T = 0.5_pr*((1._pr+a)*o12+a*b2*o22+2._pr*a*b*omega1*omega2*cos(t12))
    U = -gama*((1._pr+a)*cos(tita1)+a*b*cos(tita2))



end subroutine energy_dble_pendulum


subroutine Hamiltonian_PE(p1,p2,q1,q2,a,dp1,dp2,dq1,dq2)
    implicit none
    real(pr), intent(in)  :: a,p1,p2,q1,q2
    real(pr), intent(out) :: dp1,dp2,dq1,dq2
    
    dp1 = -q1*(1._pr+2._pr*a*q2*q2)
    dp2 = -q2*(1._pr+2._pr*a*q1*q1)
    dq1 = p1
    dq2 = p2

end subroutine Hamiltonian_PE


subroutine HamiltonianE(p1,q1,p2,q2,a,E)
    implicit none
    real(pr), intent(in)    :: a,p1,q1,p2,q2
    real(pr), intent(out)   :: E
    real(pr)              :: p12,p22,q12,q22
    
    p12 = p1*p1
    p22 = p2*p2
    q12 = q1*q1
    q22 = q2*q2

    E = 0.5_pr*(p12+p22)+0.5_pr*(q12+q22)+a*q12*q22

end subroutine HamiltonianE


end module mechanics

