module integrate
    use precision, only: pr => dp
    use mechanics
    implicit none 

    contains

    subroutine rk4_double_pendulum(ti,tf,x1_0,x2_0,x3_0,x4_0,n,tita1,tita2,omega1,omega2,a,b,gama)
        implicit none
        integer                             :: n,i
        real(pr), intent(in)                :: x1_0,x2_0,x3_0,x4_0,a,b,gama
        real(pr), allocatable, intent(out)  :: tita1(:),tita2(:),omega1(:),omega2(:)
        real(pr), intent(in)                :: ti,tf
        real(pr)                            :: k11,k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44
        real(pr)                            :: h,w1,w2,dw1,dw2

        h = (tf-ti)*(1._pr/real(n,pr))
        allocate(tita1(n+1),tita2(n+1),omega1(n+1),omega2(n+1))

        tita1(1) = x1_0
        tita2(1) = x2_0
        omega1(1) = x3_0
        omega2(1) = x4_0


        do i = 1,n
            call Lagrangian(a,b,gama,tita1(i),tita2(i),omega1(i),omega2(i),w1,w2,dw1,dw2)

            k11 = w1
            k12 = w2
            k13 = dw1
            k14 = dw2

            tita1(i+1) = tita1(i) + k11*h*0.5_pr
            tita2(i+1) = tita2(i) + k12*h*0.5_pr    
            omega1(i+1) = omega1(i) + k13*h*0.5_pr
            omega2(i+1) = omega2(i) + k14*h*0.5_pr

            call Lagrangian(a,b,gama,tita1(i+1),tita2(i+1),omega1(i+1),omega2(i+1),w1,w2,dw1,dw2)

            k21 = w1
            k22 = w2
            k23 = dw1
            k24 = dw2

            tita1(i+1) = tita1(i) + k21*h*0.5_pr
            tita2(i+1) = tita2(i) + k22*h*0.5_pr
            omega1(i+1) = omega1(i) + k23*h*0.5_pr
            omega2(i+1) = omega2(i) + k24*h*0.5_pr

            call Lagrangian(a,b,gama,tita1(i+1),tita2(i+1),omega1(i+1),omega2(i+1),w1,w2,dw1,dw2)

            k31 = w1
            k32 = w2
            k33 = dw1
            k34 = dw2

            tita1(i+1) = tita1(i) + k31*h
            tita2(i+1) = tita2(i) + k32*h
            omega1(i+1) = omega1(i) + k33*h
            omega2(i+1) = omega2(i) + k34*h

            call Lagrangian(a,b,gama,tita1(i+1),tita2(i+1),omega1(i+1),omega2(i+1),w1,w2,dw1,dw2)

            k41 = w1
            k42 = w2
            k43 = dw1
            k44 = dw2

            tita1(i+1) = tita1(i) + (k11+2._pr*(k21+k31)+k41)*h*(1._pr/6._pr)
            tita2(i+1) = tita2(i) + (k12+2._pr*(k22+k32)+k42)*h*(1._pr/6._pr)
            omega1(i+1) = omega1(i) + (k13+2._pr*(k23+k33)+k43)*h*(1._pr/6._pr)
            omega2(i+1) = omega2(i) + (k14+2._pr*(k24+k34)+k44)*h*(1._pr/6._pr)


        enddo

    end subroutine rk4_double_pendulum

    subroutine rk4_tflip(ti,tf,x1_0,x2_0,x3_0,x4_0,n,a,b,gama,tflip)
        implicit none
        integer, parameter                  :: pi=4._pr*atan(1._pr)
        integer                             :: n,i
        real(pr), intent(out)               :: tflip
        real(pr), intent(in)                :: x1_0,x2_0,x3_0,x4_0,a,b,gama
        real(pr), allocatable               :: tita1(:),tita2(:),omega1(:),omega2(:),t(:)
        real(pr), intent(in)                :: ti,tf
        real(pr)                            :: k11,k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44
        real(pr)                            :: h,w1,w2,dw1,dw2
        

        h = (tf-ti)*(1._pr/real(n,pr))
        allocate(tita1(n+1),tita2(n+1),omega1(n+1),omega2(n+1))

        tita1(1) = x1_0
        tita2(1) = x2_0
        omega1(1) = x3_0
        omega2(1) = x4_0


        t = ti + h*(/(real(i,pr), i = 0,n)/)


        do1: do i = 1,n
            if (int(tita1(i)/pi) > 1 .or. int(tita2(i)/pi) > 1) then
                tflip = t(i)
                exit do1
            else if(i == n) then
                tflip = tf + 1
            else
                call Lagrangian(a,b,gama,tita1(i),tita2(i),omega1(i),omega2(i),w1,w2,dw1,dw2)

                k11 = w1
                k12 = w2
                k13 = dw1
                k14 = dw2

                tita1(i+1) = tita1(i) + k11*h*0.5_pr
                tita2(i+1) = tita2(i) + k12*h*0.5_pr    
                omega1(i+1) = omega1(i) + k13*h*0.5_pr
                omega2(i+1) = omega2(i) + k14*h*0.5_pr

                call Lagrangian(a,b,gama,tita1(i+1),tita2(i+1),omega1(i+1),omega2(i+1),w1,w2,dw1,dw2)

                k21 = w1
                k22 = w2
                k23 = dw1
                k24 = dw2

                tita1(i+1) = tita1(i) + k21*h*0.5_pr
                tita2(i+1) = tita2(i) + k22*h*0.5_pr
                omega1(i+1) = omega1(i) + k23*h*0.5_pr
                omega2(i+1) = omega2(i) + k24*h*0.5_pr

                call Lagrangian(a,b,gama,tita1(i+1),tita2(i+1),omega1(i+1),omega2(i+1),w1,w2,dw1,dw2)

                k31 = w1
                k32 = w2
                k33 = dw1
                k34 = dw2

                tita1(i+1) = tita1(i) + k31*h
                tita2(i+1) = tita2(i) + k32*h
                omega1(i+1) = omega1(i) + k33*h
                omega2(i+1) = omega2(i) + k34*h

                call Lagrangian(a,b,gama,tita1(i+1),tita2(i+1),omega1(i+1),omega2(i+1),w1,w2,dw1,dw2)

                k41 = w1
                k42 = w2
                k43 = dw1
                k44 = dw2

                tita1(i+1) = tita1(i) + (k11+2._pr*(k21+k31)+k41)*h*(1._pr/6._pr)
                tita2(i+1) = tita2(i) + (k12+2._pr*(k22+k32)+k42)*h*(1._pr/6._pr)
                omega1(i+1) = omega1(i) + (k13+2._pr*(k23+k33)+k43)*h*(1._pr/6._pr)
                omega2(i+1) = omega2(i) + (k14+2._pr*(k24+k34)+k44)*h*(1._pr/6._pr)
            endif
        enddo do1

    end subroutine rk4_tflip








    subroutine rk4_Hamiltoniano_PE(ti,tf,p1_0,p2_0,q1_0,q2_0,n,p1,p2,q1,q2,a)
        implicit none
        integer                             :: n,i
        real(pr), intent(in)                :: p1_0,p2_0,q1_0,q2_0,a
        real(pr), allocatable, intent(out)  :: p1(:),p2(:),q1(:),q2(:)
        real(pr), intent(in)                :: ti,tf
        real(pr)                            :: k11,k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44
        real(pr)                            :: h,dp1,dp2,dq1,dq2

        h = (tf-ti)*(1._pr/real(n,pr))
        allocate(p1(n+1),p2(n+1),q1(n+1),q2(n+1))

        p1(1) = p1_0
        p2(1) = p2_0
        q1(1) = q1_0
        q2(1) = q2_0


        do i = 1,n
            call Hamiltonian_PE(p1(i),p2(i),q1(i),q2(i),a,dp1,dp2,dq1,dq2)

            k11 = dp1
            k12 = dp2
            k13 = dq1
            k14 = dq2

            p1(i+1) = p1(i) + k11*h*0.5_pr
            p2(i+1) = p2(i) + k12*h*0.5_pr    
            q1(i+1) = q1(i) + k13*h*0.5_pr
            q2(i+1) = q2(i) + k14*h*0.5_pr


            call Hamiltonian_PE(p1(i+1),p2(i+1),q1(i+1),q2(i+1),a,dp1,dp2,dq1,dq2)

            k21 = dp1
            k22 = dp2
            k23 = dq1
            k24 = dq2

            p1(i+1) = p1(i) + k21*h*0.5_pr
            p2(i+1) = p2(i) + k22*h*0.5_pr    
            q1(i+1) = q1(i) + k23*h*0.5_pr
            q2(i+1) = q2(i) + k24*h*0.5_pr

            call Hamiltonian_PE(p1(i+1),p2(i+1),q1(i+1),q2(i+1),a,dp1,dp2,dq1,dq2)

            k31 = dp1
            k32 = dp2
            k33 = dq1
            k34 = dq2

            p1(i+1) = p1(i) + k31*h
            p2(i+1) = p2(i) + k32*h    
            q1(i+1) = q1(i) + k33*h
            q2(i+1) = q2(i) + k34*h

            call Hamiltonian_PE(p1(i+1),p2(i+1),q1(i+1),q2(i+1),a,dp1,dp2,dq1,dq2)

            k41 = dp1
            k42 = dp2
            k43 = dq1
            k44 = dq2

            p1(i+1) = p1(i) + (k11+2._pr*(k21+k31)+k41)*h*(1._pr/6._pr)
            p2(i+1) = p2(i) + (k12+2._pr*(k22+k32)+k42)*h*(1._pr/6._pr)
            q1(i+1) = q1(i) + (k13+2._pr*(k23+k33)+k43)*h*(1._pr/6._pr)
            q2(i+1) = q2(i) + (k14+2._pr*(k24+k34)+k44)*h*(1._pr/6._pr)


        enddo

    end subroutine rk4_Hamiltoniano_PE


    subroutine rk4_PEH_Poincare(ti,tf,n,np,p1_0,q1_0,p2_0,q2_0,a,p1_map,q1_map)
        implicit none
        integer                             :: n,i,np,counter
        real(pr), intent(in)                :: p1_0,p2_0,q1_0,q2_0,a
        real(pr), allocatable, intent(out)  :: p1_map(:),q1_map(:) 
        real(pr), allocatable               :: p1(:),q1(:),p2(:),q2(:)
        real(pr), intent(in)                :: ti,tf
        real(pr)                            :: k11,k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44
        real(pr)                            :: h,dp1,dp2,dq1,dq2


        h = (tf-ti)*(1._pr/real(n,pr))
        allocate(p1(n+1),p2(n+1),q1(n+1),q2(n+1))

        allocate(p1_map(np),q1_map(np))

        p1(1) = p1_0
        p2(1) = p2_0
        q1(1) = q1_0
        q2(1) = q2_0


        counter = 1

        do1: do i = 1,n
            call Hamiltonian_PE(p1(i),p2(i),q1(i),q2(i),a,dp1,dp2,dq1,dq2)

            k11 = dp1
            k12 = dp2
            k13 = dq1
            k14 = dq2

            p1(i+1) = p1(i) + k11*h*0.5_pr
            p2(i+1) = p2(i) + k12*h*0.5_pr    
            q1(i+1) = q1(i) + k13*h*0.5_pr
            q2(i+1) = q2(i) + k14*h*0.5_pr


            call Hamiltonian_PE(p1(i+1),p2(i+1),q1(i+1),q2(i+1),a,dp1,dp2,dq1,dq2)

            k21 = dp1
            k22 = dp2
            k23 = dq1
            k24 = dq2

            p1(i+1) = p1(i) + k21*h*0.5_pr
            p2(i+1) = p2(i) + k22*h*0.5_pr    
            q1(i+1) = q1(i) + k23*h*0.5_pr
            q2(i+1) = q2(i) + k24*h*0.5_pr

            call Hamiltonian_PE(p1(i+1),p2(i+1),q1(i+1),q2(i+1),a,dp1,dp2,dq1,dq2)

            k31 = dp1
            k32 = dp2
            k33 = dq1
            k34 = dq2

            p1(i+1) = p1(i) + k31*h
            p2(i+1) = p2(i) + k32*h    
            q1(i+1) = q1(i) + k33*h
            q2(i+1) = q2(i) + k34*h

            call Hamiltonian_PE(p1(i+1),p2(i+1),q1(i+1),q2(i+1),a,dp1,dp2,dq1,dq2)

            k41 = dp1
            k42 = dp2
            k43 = dq1
            k44 = dq2

            p1(i+1) = p1(i) + (k11+2._pr*(k21+k31)+k41)*h*(1._pr/6._pr)
            p2(i+1) = p2(i) + (k12+2._pr*(k22+k32)+k42)*h*(1._pr/6._pr)
            q1(i+1) = q1(i) + (k13+2._pr*(k23+k33)+k43)*h*(1._pr/6._pr)
            q2(i+1) = q2(i) + (k14+2._pr*(k24+k34)+k44)*h*(1._pr/6._pr)

            if (i > int(n*0.5_pr)) then
                if (q2(i)*q2(i+1) < 0._pr) then
                    q1_map(counter) = (q1(i)+q1(i+1))*(1._pr/2._pr)
                    p1_map(counter) = (p1(i)+p1(i+1))*(1._pr/2._pr)
                    counter = counter + 1
                    if (counter == np + 1) then
                        exit do1
                    endif
                endif 
            endif
        enddo do1

    end subroutine rk4_PEH_Poincare




end module integrate
