!------------------------------------------------------------------------------
! Institution, Affiliation
!------------------------------------------------------------------------------
!
! MODULE:  Module name
!
!> @author
!> Author Name}
!
! DESCRIPTION: 
!>  Short module description
!
! REVISION HISTORY:
! dd Mmm yyyy - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module module_EDO_segundo_orden_flip
    use module_presition
    use module_double_pendulum

    implicit none
    
contains
    subroutine RK4_four_eq_flip(RK4_type,n,a,b,y1_0,y2_0,y3_0,y4_0,x_flip,function_type,input_type)
        implicit none

        ! Data dictionary: declare calling parameter types & definitions
        integer(sp),            intent(in) :: RK4_type
        integer(sp),            intent(in) :: n              ! points
        integer(sp),            intent(in) :: function_type  ! tipo de función a integrar
        integer(sp),            intent(in) :: input_type     ! type of input to create function
        real(dp),               intent(in) :: y1_0           ! initial condition
        real(dp),               intent(in) :: y2_0           ! initial condition
        real(dp),               intent(in) :: y3_0           ! initial condition
        real(dp),               intent(in) :: y4_0           ! initial condition
        real(dp),               intent(in) :: a,b            ! intervalo de validez para la solución aproximada
        real(dp),               intent(out) :: x_flip        ! flip coordinate
        
        ! Data dictionary: declare local variables types & definitions
        integer(sp) 			:: i,index_prev 			! index loop
        real(dp) 				:: k1_1, k1_2, k1_3, k1_4	! recurrence relations
        real(dp) 				:: k2_1, k2_2, k2_3, k2_4	! recurrence relations
        real(dp) 				:: k3_1, k3_2, k3_3, k3_4	! recurrence relations
        real(dp) 				:: k4_1, k4_2, k4_3, k4_4	! recurrence relations
        real(dp) 				:: y1_improved_k2			! improved function evaluation
        real(dp) 				:: y1_improved_k3			! improved function evaluation
        real(dp) 				:: y1_improved_k4			! improved function evaluation
        real(dp) 				:: y2_improved_k2			! improved function evaluation
        real(dp) 				:: y2_improved_k3			! improved function evaluation
        real(dp) 				:: y2_improved_k4			! improved function evaluation
        real(dp) 				:: y3_improved_k2			! improved function evaluation
        real(dp) 				:: y3_improved_k3			! improved function evaluation
        real(dp) 				:: y3_improved_k4			! improved function evaluation
        real(dp) 				:: y4_improved_k2			! improved function evaluation
        real(dp) 				:: y4_improved_k3			! improved function evaluation
        real(dp) 				:: y4_improved_k4			! improved function evaluation
        real(dp) 				:: h_improved 				! improved step evaluation
        real(dp), dimension(n) 	:: x 						! grid points vector
        real(dp) 				:: h 						! step
        real(dp) 				:: func_lag1,func_lag2,func_lag3,func_lag4 				! custom functions
        real(dp), parameter     :: pi=4._dp*atan(1._dp)
        real(dp), parameter     :: a1=1._dp,b1=1._dp,g1=1._dp
        real(dp), allocatable   :: y1_RK4(:), y2_RK4(:)
        real(dp), allocatable   :: y3_RK4(:), y4_RK4(:)


        allocate(y1_RK4(n),y2_RK4(n),y3_RK4(n),y4_RK4(n))        
        select case(input_type)
            case(1) ! using lagrangian_dble_pendulum(q1,q2,dq1,dq2,a,b,c,w1,w2,dw1,dw2,function_type)
                select case (RK4_type)
                    case(1) ! Método clásico

                        ! remember --> a=alpha=m2/m1; b=beta=l2/l1; c=gamma=g/l1

                        h = abs(b - a) * ( 1._dp / (real(n,dp)-1._dp) )
                        x(1) = a


                        y1_RK4(1) = y1_0 ! initial condition of generalized coordinate 1
                        y2_RK4(1) = y2_0 ! initial condition of generalized coordinate 2
                        y3_RK4(1) = y3_0 ! initial condition of generalized velocity 1
                        y4_RK4(1) = y4_0 ! initial condition of generalized velocity 2


                        do i = 2, n, 1
                            if (y1_RK4(i) <= pi .or. y2_RK4(i) <= pi ) then
                                index_prev = i-1
                                x(i) = x(index_prev) + h
                        
                                call lagrangian_dble_pendulum(y1_RK4(index_prev),y2_RK4(index_prev),y3_RK4(index_prev),&
                                y4_RK4(index_prev),a1,b1,g1,func_lag1,func_lag2,func_lag3,func_lag4,1_sp)
                        
                                k1_1 = func_lag1
                                k2_1 = func_lag2
                                k3_1 = func_lag3
                                k4_1 = func_lag4
                                
                                h_improved = 0.5_dp*h
                                
                                y1_improved_k2 = y1_RK4(index_prev) + (k1_1*h_improved)
                                y2_improved_k2 = y2_RK4(index_prev) + (k2_1*h_improved)
                                y3_improved_k2 = y3_RK4(index_prev) + (k3_1*h_improved)
                                y4_improved_k2 = y4_RK4(index_prev) + (k4_1*h_improved)
                                
                                call lagrangian_dble_pendulum(y1_improved_k2,y2_improved_k2,y3_improved_k2,&
                                y4_improved_k2,a1,b1,g1,func_lag1,func_lag2,func_lag3,func_lag4,1_sp)
                                
                                k1_2 = func_lag1
                                k2_2 = func_lag2
                                k3_2 = func_lag3
                                k4_2 = func_lag4
                                
                                y1_improved_k3 = y1_RK4(index_prev) + (k1_2*h_improved)
                                y2_improved_k3 = y2_RK4(index_prev) + (k2_2*h_improved)
                                y3_improved_k3 = y3_RK4(index_prev) + (k3_2*h_improved)
                                y4_improved_k3 = y4_RK4(index_prev) + (k4_2*h_improved)
                                
                                call lagrangian_dble_pendulum(y1_improved_k3,y2_improved_k3,y3_improved_k3,&
                                y4_improved_k3,a1,b1,g1,func_lag1,func_lag2,func_lag3,func_lag4,1_sp)
                                k1_3 = func_lag1
                                k2_3 = func_lag2
                                k3_3 = func_lag3
                                k4_3 = func_lag4
                                
                                y1_improved_k4 = y1_RK4(index_prev) + (k1_3*h)
                                y2_improved_k4 = y2_RK4(index_prev) + (k2_3*h)
                                y3_improved_k4 = y3_RK4(index_prev) + (k3_3*h)
                                y4_improved_k4 = y4_RK4(index_prev) + (k4_3*h)
                                
                                call lagrangian_dble_pendulum(y1_improved_k4,y2_improved_k4,y3_improved_k4,&
                                y4_improved_k4,a1,b1,g1,func_lag1,func_lag2,func_lag3,func_lag4,1_sp)
                                k1_4 = func_lag1
                                k2_4 = func_lag2
                                k3_4 = func_lag3
                                k4_4 = func_lag4
                                
                                ! y(i+1) = (y(i) + increment_function)
                                y1_RK4(i) = y1_RK4(index_prev) + (1._dp/6._dp)*h*(k1_1 + 2._dp*(k1_2+k1_3)+k1_4)
                                y2_RK4(i) = y2_RK4(index_prev) + (1._dp/6._dp)*h*(k2_1 + 2._dp*(k2_2+k2_3)+k2_4)
                                y3_RK4(i) = y3_RK4(index_prev) + (1._dp/6._dp)*h*(k3_1 + 2._dp*(k3_2+k3_3)+k3_4)
                                y4_RK4(i) = y4_RK4(index_prev) + (1._dp/6._dp)*h*(k4_1 + 2._dp*(k4_2+k4_3)+k4_4)
                            else
                                x_flip = x(i) ! detected flip
                                write(*,*) i, x_flip, y1_0, y2_0, y1_RK4(i), y2_RK4(i)
                                exit
                            end if

                            x_flip = (b + 1._dp) ! no-detected flip
                        end do
                        
                    case default
                        write(*,*) 'Invalid RK4 type'
                end select
            !case(2) ! add to use f_2D(x,y,function_type)
            !case(3) ! add to use f_1D(x,function_type)
            case default
                write(*,*) 'Invalid input type'
        end select
    deallocate(y1_RK4,y2_RK4,y3_RK4,y4_RK4) 
    end subroutine RK4_four_eq_flip
end module module_EDO_segundo_orden_flip