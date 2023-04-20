module miscelaneous
    use precision
    implicit none

    contains

    subroutine exp_fit(x,y,m,b)
        use precision, only: dp=>dp
        implicit none

        real(dp), intent(in)  :: y(:),x(:)
        real(dp), allocatable :: v(:),u(:)
        real(dp), intent(out) :: m,b
        real(dp)              :: v_m,u_m,u2_m,uv_m
        integer(dp)           :: i,n

        n = size(y)

        allocate(v(n),u(n))

        u_m=0._dp;v_m=0._dp;uv_m=0._dp;u2_m=0._dp

        v = log(y);u = log(x)

        u_m = sum(u)/real(n,dp)
        v_m = sum(v)/real(n,dp)

        do i = 1,n
            uv_m = uv_m + u(i)*v(i)
            u2_m = u2_m + u(i)*u(i)
        enddo

        uv_m = uv_m/real(n,dp)
        u2_m = u2_m/real(n,dp)

        m = (uv_m - u_m*v_m)/(u2_m - u_m*u_m)
        b = exp(v_m - m*u_m)

    end subroutine exp_fit

    function nsphere(x,d)
        use precision, only: dp=>dp,prs=>sp
        implicit none
        integer(dp)  :: d
        real(dp)     :: x,nsphere

        nsphere = sqrt((1._dp - x*x)**real(d-1,dp))
    end function nsphere

    function factorial(n)
        use precision, only: dp=>dp,prs=>sp
        implicit none
        real(dp)    :: n,factorial
        integer(dp) :: i,n_loop

        n_loop = int(n,dp)
        factorial = 1_dp
        do i = n_loop,1,-1
            factorial = factorial*real(i,dp)
        enddo
    end function factorial

    function vol_sphere(d)
        use precision, only: dp=>dp,prs=>sp
        implicit none
        real(dp), parameter  :: pi=4._dp*atan(1._dp)
        real(dp)             :: vol_sphere,d

        if (mod(int(d,dp),2) == 0) then
            vol_sphere = pi**(d/2._dp)*(1._dp/factorial(d*0.5_dp))
        else
            vol_sphere = (2_dp**d)*pi**((d-1_dp)/2_dp)*factorial((d-1_dp)/2_dp)*(1._dp/factorial(d))
        endif
    end function vol_sphere

    function gauss_norm(x,sigma)
        use precision
        implicit none

        real(dp)            :: x,sigma
        real(dp), parameter :: pi=4._dp*atan(1._dp)
        real(dp)            :: gauss_norm,factor_1,factor_2

        factor_1 = 1._dp/(sqrt(2._dp*pi)*sigma)
        factor_2 = x*(1._dp/sigma)
        gauss_norm = factor_1*exp(-0.5_dp*factor_2*factor_2)
    end function gauss_norm

    subroutine mc_step(T,U,states,seed)
        use precision,prs=>sp
        use mtmod ,only: sgrnd,grnd
        implicit none

        integer(dp), intent(inout) :: states(:,:)
        real(dp), intent(in)       :: T
        real(dp), intent(out)      :: U
        integer(4),intent(inout)   :: seed!,seed_n(8)
        integer(dp)                :: i,j,k,nn
        real(dp)                   :: vecinos,delta_U,Bolt,ran_u

        !call date_and_time(values=seed_n)
        !seed = 9284!seed_n(8)*seed_n(7)*seed_n(6)+seed_n(5)
        call sgrnd(seed)

        nn = size(states(:,1))-2_dp

        U = 0._dp

        do k = 1,nn*nn
            ! Generamos coordenadas aleatorias en nn x nn
            ran_u = grnd()*real(nn,dp)
            i = int(ran_u)+2_dp
            ran_u = grnd()*real(nn,dp)
            j = int(ran_u)+2_dp
            if (j == nn+2_dp) then
                j = nn+1_dp
            else if (i == nn+2_dp) then
                i = nn+1_dp
            endif
            ! Calculamos los primeros vecinos
            vecinos = real(states(i-1,j)+states(i+1,j)+&
            & states(i,j-1)+states(i,j+1),dp)
            ! Calculamos el delta de energía
            delta_U = 2._dp*real(states(i,j),dp)*vecinos
            ! Aceptamos flip si delta es menor o igual a cero.
            if (delta_U <= 0._dp) then
                states(i,j) = -1_dp*states(i,j)
                ! De ser necesario actualizamos condiciones de borde
                if (i == 2_dp .or. j == 2_dp .or. i == nn+1_dp .or. j == nn+1_dp) then
                    states(1,j) = states(nn+1,j)
                    states(nn+2,j) = states(2,j)
                    states(i,1) = states(i,nn+1)
                    states(i,nn+2) = states(i,2)
                endif
                U = U + delta_U
            ! En caso de ser positivo aplicamos probabilidad de Boltzmann
            else
                ran_u = grnd()
                ! Si T == 0, la exponencial
                if (T == 0._dp) then
                Bolt = 0._dp
                else
                    Bolt = exp(-delta_U/T)
                endif
                ! Si la distribución de Boltzmann arroja el número mayor, aceptamos el flip
                if (Bolt >= ran_u) then
                    states(i,j) = -1_dp*states(i,j)
                    ! Si es necesario actualizamos condiciones de borde
                    if (i == 2_dp .or. j == 2_dp .or. i == nn+1_dp .or. j == nn+1_dp) then
                        states(1,j) = states(nn+1,j)
                        states(nn+2,j) = states(2,j)
                        states(i,1) = states(i,nn+1)
                        states(i,nn+2) = states(i,2)
                    endif
                    U = U + delta_U
                else
                    U = U
                endif
            endif
        enddo
        seed = abs(int(grnd()*9000))+150


    end subroutine mc_step
   

    subroutine fill_random(states)
        use precision,prs=>sp
        use mtmod
        implicit none

        integer(dp), intent(inout) :: states(:,:)
        integer(dp)                :: i,j,nn
        integer(prs)                :: seed,seed_n(8)
        real(dp)                   :: ran_u

        nn = size(states(1,:))-2_dp

        call date_and_time(values=seed_n)
        seed = seed_n(8)*seed_n(7)*seed_n(6)+seed_n(5)
        call sgrnd(seed)

        ! LLenamos matriz de forma aleatoria
        do i = 2,nn+1
            do j = 2,nn+1
                ran_u = grnd()
                if (ran_u <= 0.5_dp) then
                    states(i,j) = 1_dp
                else
                    states(i,j) = -1_dp
                endif
            enddo
        enddo
        ! Aplicamos condiciones de borde periódicas
        states(1,:) = states(nn+1,:)
        states(nn+2,:) = states(2,:)
        states(:,1) = states(:,nn+1)
        states(:,nn+2) = states(:,2)

    end subroutine fill_random


    function H(states)
        use precision
        implicit none

        integer(dp)     :: ii,jj,n_states
        real(dp)        :: H,vec
        integer(dp)     :: states(:,:)

        n_states = size(states(1,:))-2_dp

        H = 0._dp

        do ii = 2,n_states + 1
            do jj = 2,n_states + 1
                vec = real(states(ii-1,jj)+states(ii+1,jj)+&
                & states(ii,jj-1)+states(ii,jj+1),dp)
                H = H + - (real(states(ii,jj),dp)*vec)
            enddo
        enddo

        H = H*0.5_dp

    end function H
    
    function M(states)
        use precision,sp=>sp
        implicit none
    
        integer(dp)     :: ii,jj,n_states
        real(dp)        :: M
        integer(dp)     :: states(:,:)
    
        n_states = size(states(1,:))-2_dp
    
        M = 0._dp
    
        do ii = 2,n_states + 1
            do jj = 2,n_states + 1
                M = M + states(ii,jj)
            enddo
        enddo

        M = M*(1._dp/real(n_states*n_states,dp))

        M = abs(M)
    
    end function M

    function z(T_ad)
        use precision
        implicit none

        real(dp), parameter    :: kb=1.38066e-23
        real(dp)               :: z,T_ad

        z = exp(-2._dp/T_ad)
    end function z

    function M_exact(T_ad)
        use precision
        implicit none

        real(dp), parameter  :: kb=1.38066e-23,Tc = 2.2676_dp
        real(dp)             :: M_exact,zz,T_ad

        if (T_ad > Tc) then
            M_exact = 0._dp
        else
            zz = z(T_ad)
            M_exact = (((1+zz*zz)**0.25_dp)*&
            &(1._dp-6._dp*zz*zz+zz*zz*zz*zz)**0.125_dp)/(sqrt(1-zz*zz))
        endif
    end function M_exact
    
    function media(x)
        use precision
        implicit none

        real(dp)       :: x(:)
        real(dp)       :: media,nn

        nn = real(size(x),dp)

        media = sum(x)/nn

    end function media

    function varianza(x)
        use precision
        implicit none

        real(dp)       :: x(:)
        real(dp)       :: varianza,media_x
        integer(dp)    :: i,nn

        nn = size(x)

        media_x = sum(x)/(real(nn,dp))
        varianza = 0._dp

        do i = 1,nn
            varianza = varianza + (x(i)-media_x)*(x(i)-media_x)
        enddo

        varianza = varianza/(real(nn-1,dp))

    end function varianza

    subroutine thermalization(T,MCS,trans,ne)
        use precision
        implicit none

        real(dp), intent(inout)   :: T
        integer(dp), intent(out)  :: MCS,trans,ne

        if (T >= 1.75_dp .and. T < 3.45_dp) then
            MCS = 20000_dp
            trans = 10000_dp
            ne = 10_dp
            T = T + 0.02
        else if (T < 1.75_dp) then
            MCS = 5000_dp
            trans = 3000_dp
            ne = 10_dp
            T = T + 0.2_dp
        else
            MCS = 10000_dp
            trans = 5000_dp
            ne = 10_dp
            T = T + 0.2_dp
        endif
    
    end subroutine thermalization


    subroutine autocorrelation(i,T_corr,U_i,M_i,U_v,M_v,t_U,t_M,u_shift,m_shift,n_points)
        use precision, only: dp=>dp
        implicit none
    
        integer(dp), intent(in)    :: i,T_corr
        real(dp), intent(in)       :: U_i,M_i
        real(dp), intent(inout)    :: U_v(:),M_v(:),t_U(:),t_M(:),u_shift(:),m_shift(:)
        integer(dp), intent(inout)    :: n_points(:)
        integer(dp)                :: n,j
        integer(dp), save          :: corr_idx
    
        n = size(U_v)

        if (i == T_corr + 1_dp) then
            corr_idx = 1_dp
            u_shift = cshift(U_v,1);m_shift = cshift(M_v,1)
            u_shift(n) = U_i;m_shift(n) = M_i
        else if (mod(i-1_dp,T_corr) == 0_dp .and. i > T_corr + 1_dp) then  
            corr_idx = 1_dp
            U_v = u_shift;M_v=m_shift
            u_shift = cshift(U_v,1);m_shift = cshift(M_v,1)
            u_shift(n) = U_i;m_shift(n) = M_i
        else
            corr_idx = corr_idx + 1_dp
            u_shift = cshift(u_shift,1);m_shift = cshift(m_shift,1)
            u_shift(n) = U_i;m_shift(n) = M_i
        endif
        
        if (i <= T_corr*1000_dp) then
            n_points(corr_idx) = n_points(corr_idx) + T_corr
        else if (i > T_corr*1000_dp) then
            n_points(corr_idx) = n_points(corr_idx) + T_corr - corr_idx
        endif
        t_U(corr_idx) = t_U(corr_idx) + dot_product(U_v,u_shift)
        t_M(corr_idx) = t_M(corr_idx) + dot_product(M_v,m_shift)
        
        

    end subroutine autocorrelation
    
    subroutine sc_fill(r,N,a)
        use precision
        implicit none

        integer(dp)                         :: i,j,k,L,idx
        real(dp), intent(in)                :: N,a
        real(dp), allocatable, intent(out)  :: r(:,:)

        L = anint((N)**(1._dp/3._dp))

        allocate(r(int(N),3_dp))

        do i = 1,L
            do j = 1,L
                do k = 1,L
                    idx = (i-1)*L*L + (j-1)*L + k
                    r(idx,1) = a*real(i,dp)+0.5_dp*a
                    r(idx,2) = a*real(j,dp)+0.5_dp*a
                    r(idx,3) = a*real(k,dp)+0.5_dp*a
                enddo
            enddo
        enddo

    end subroutine sc_fill

    subroutine fcc_fill(r,N,rho,L)
        use precision
        implicit none

        integer(dp)             :: i,j,k,idx,n_cells
        real(dp), intent(in)    :: rho
        integer(dp), intent(in) :: N
        real(dp)                :: L_2
        real(dp), intent(out)   :: L
        real(dp), allocatable   :: r(:,:)

        n_cells = int(anint((real(N,dp)*0.25_dp)**(1._dp/3._dp)),dp)
        L = (4._dp/rho)**(1._dp/3._dp)
        L_2 = L*0.5_dp

        allocate(r(N,3_dp))
        r = 0._dp

        idx = 0_dp  

        do i = 0,n_cells-1
            do j = 0,n_cells-1
                do k = 0,n_cells-1
                    idx = idx + 1_dp
                    r(idx,:) = [L*real(i,dp),L*real(j,dp),L*real(k,dp)]
                    idx = idx + 1_dp
                    r(idx,:) = [L*real(i,dp)+L_2,L*real(j,dp)+L_2,L*real(k,dp)]
                    idx = idx + 1_dp
                    r(idx,:) = [L*real(i,dp),L*real(j,dp)+L_2,L*real(k,dp)+L_2]
                    idx = idx + 1_dp
                    r(idx,:) = [L*real(i,dp)+L_2,L*real(j,dp),L*real(k,dp)+L_2]
                enddo
            enddo
        enddo
        L = real(n_cells,dp)*L
        
        r = r - L*0.5_dp
    end subroutine fcc_fill

    subroutine md_pbc(r,L)
        use precision
        implicit none

        integer(dp)             :: i
        real(dp), intent(inout) :: r(:,:)
        real(dp), intent(in)    :: L

        do i = 1,size(r(:,1))
            r(i,1) = r(i,1) - L*anint(r(i,1)/L)
            r(i,2) = r(i,2) - L*anint(r(i,2)/L)
            r(i,3) = r(i,3) - L*anint(r(i,3)/L)
        enddo


    end subroutine md_pbc

    subroutine init_gaussian_v(v,N,seed)
        use precision
        use rangens
        implicit none

        integer(dp)                         :: i
        integer(dp), intent(in)             :: N
        integer(sp), intent(in)             :: seed
        real(dp), intent(inout)             :: v(N,3)
        real(dp)                            :: v1,v2,v3

        v = 0._dp

        do i = 1,N
            v(i,1) = gaussdev(seed)
            v(i,2) = gaussdev(seed)
            v(i,3) = gaussdev(seed)
        enddo

        v1 = sum(v(:,1))/real(N,dp)
        v2 = sum(v(:,2))/real(N,dp)
        v3 = sum(v(:,3))/real(N,dp)

        v(:,1) = v(:,1) - v1
        v(:,2) = v(:,2) - v2
        v(:,3) = v(:,3) - v3

    end subroutine init_gaussian_v


    subroutine init_uniform_v(v,N,seed)
        use precision
        use mtmod, only: sgrnd,grnd
        implicit none

        integer(dp)                         :: i
        integer(dp), intent(in)             :: N
        integer(sp), intent(in)             :: seed
        real(dp), intent(inout)             :: v(N,3)
        real(dp)                            :: v1,v2,v3

        v = 0._dp

        call sgrnd(seed)

        do i = 1,N
            v(i,1) = -0.5_dp + grnd()
            v(i,2) = -0.5_dp + grnd()
            v(i,3) = -0.5_dp + grnd()
        enddo

        v1 = sum(v(:,1))/real(N,dp)
        v2 = sum(v(:,2))/real(N,dp)
        v3 = sum(v(:,3))/real(N,dp)

        v(:,1) = v(:,1) - v1
        v(:,2) = v(:,2) - v2
        v(:,3) = v(:,3) - v3

    end subroutine init_uniform_v


    subroutine U_LJ(r,U)
        use precision
        implicit none

        real(dp), intent(in)    :: r
        real(dp), intent(out)   :: U
        real(dp)                :: r_1,r6,r12

        r_1 = (1._dp/r)
        r6 = r_1**6

        r12 = r_1**12

        U = 4._dp*(r12-r6)

    end subroutine U_LJ

    subroutine U_ij(L,jdx,idx,rcut,x,Uscale,U)
        use precision
        implicit none

        integer(dp), intent(in) :: jdx,idx
        real(dp), intent(in)    :: x(:,:),L,rcut,Uscale
        real(dp), intent(out)   :: U
        real(dp)                :: Uij,dx1,dx2,dx3,rij


        Uij = 0._dp; U = 0._dp
        dx1 = x(jdx,1)-x(idx,1)
        dx2 = x(jdx,2)-x(idx,2)
        dx3 = x(jdx,3)-x(idx,3)
        dx1 = dx1 - L*anint(dx1/L)
        dx2 = dx2 - L*anint(dx2/L)
        dx3 = dx3 - L*anint(dx3/L)
        rij = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3)
        if (rij <= rcut) then
            call U_LJ(rij,Uij)
            U = Uij + Uscale
        else
            U = 0._dp
        endif
    
    end subroutine U_ij


    subroutine f_LJ(r,u)
        use precision
        implicit none

        real(dp), intent(in)   :: r
        real(dp), intent(out)  :: u
        real(dp)               :: r6,r12

        r6 = (1._dp/r)
        r6 = r6*r6*r6*r6*r6*r6
        r12 = r6*r6
        u = 48._dp*(r12-0.5_dp*r6)/(r*r)

    end subroutine f_LJ

    subroutine f_ij(L,jdx,idx,rcut,x,f,w)
        use precision
        implicit none

        integer(dp), intent(in) :: jdx,idx
        real(dp), intent(in)    :: x(:,:),L,rcut
        real(dp), intent(inout) :: f(:,:)
        real(dp), intent(out)   :: w
        real(dp)                :: u,dx1,dx2,dx3,rij

        dx1 = x(jdx,1)-x(idx,1)
        dx2 = x(jdx,2)-x(idx,2)
        dx3 = x(jdx,3)-x(idx,3)
        dx1 = dx1 - L*anint(dx1/L)
        dx2 = dx2 - L*anint(dx2/L)
        dx3 = dx3 - L*anint(dx3/L)
        rij = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3)
        if (rij <= rcut) then
            call f_LJ(rij,u)
            w = w + rij*rij*u
            f(jdx,1) = f(jdx,1) + dx1*u
            f(jdx,2) = f(jdx,2) + dx2*u
            f(jdx,3) = f(jdx,3) + dx3*u
            f(idx,1) = f(idx,1) - dx1*u
            f(idx,2) = f(idx,2) - dx2*u
            f(idx,3) = f(idx,3) - dx3*u
        endif
        
    end subroutine f_ij

    subroutine system_forces(r,rcut,f,L,w)
        use precision
        implicit none

        real(dp), intent(in)    :: r(:,:)
        real(dp), intent(in)    :: rcut,L
        real(dp), intent(inout) :: f(:,:)
        real(dp), intent(out)   :: w
        integer(dp)             :: i,j
        real(dp)                :: rij,dx1,dx2,dx3,u

        f = 0._dp
        w = 0._dp
    
        do j = 2,size(r(:,1))
            do i = 1,j-1
                call f_ij(L,j,i,rcut,r,f,w)
            enddo
        enddo
    end subroutine system_forces

    function icell(ix,iy,iz,M_cell)
        use precision
        implicit none

        integer(dp) :: ix,iy,iz,M_cell,icell

        icell = 1_dp + mod(ix - 1_dp + M_cell,M_cell)  &
        + mod(iy - 1_dp + M_cell,M_cell)*M_cell &
        + mod(iz - 1 + M_cell,M_cell)*M_cell*M_cell

    end function icell

    subroutine maps(mapa,M_cell)
        use precision
        implicit none

        integer(dp), intent(in)     :: M_cell
        integer(dp), intent(out)    :: mapa(:)
        integer(dp)                 :: ix,iy,iz,imap

        do ix = 1,M_cell
            do iy = 1,M_cell
                do iz = 1,M_cell
                    imap = (icell(ix,iy,iz,M_cell) - 1._dp)*13_dp

                    mapa(imap+1) = icell(ix+1,iy,iz,M_cell)
                    mapa(imap+2) = icell(ix+1,iy+1,iz,M_cell)
                    mapa(imap+3) = icell(ix,iy+1,iz,M_cell)
                    mapa(imap+4) = icell(ix-1,iy+1,iz,M_cell)
                    mapa(imap+5) = icell(ix+1,iy,iz-1,M_cell)
                    mapa(imap+6) = icell(ix+1,iy+1,iz-1,M_cell)
                    mapa(imap+7) = icell(ix,iy+1,iz-1,M_cell)
                    mapa(imap+8) = icell(ix-1,iy+1,iz-1,M_cell)
                    mapa(imap+9) = icell(ix+1,iy,iz+1,M_cell)
                    mapa(imap+10) = icell(ix+1,iy+1,iz+1,M_cell)
                    mapa(imap+11) = icell(ix,iy+1,iz+1,M_cell)
                    mapa(imap+12) = icell(ix-1,iy+1,iz+1,M_cell)
                    mapa(imap+13) = icell(ix,iy,iz+1,M_cell)
                enddo
            enddo
        enddo    
    end subroutine maps

    subroutine links(M_cell,L,rcut,x,head,list)
        use precision
        implicit none

        integer(dp), intent(in)     :: M_cell
        real(dp), intent(in)        :: x(:,:),L,rcut
        integer(dp), intent(out)    :: head(:),list(:)
        integer(dp)                 :: idxcell,Ncell,N,i
        real(dp)                    :: Lc_inv

        Ncell = size(head)
        N = size(list)

        head(1:Ncell) = 0_dp

        Lc_inv = real(M_cell)/L
        
        do i = 1,N
            idxcell = 1_dp + int((x(i,1)+0.5_dp*L)*Lc_inv) + &
                int((x(i,2)+0.5_dp*L)*Lc_inv)*M_cell + &
                int((x(i,3)+0.5_dp*L)*Lc_inv)*M_cell*M_cell
            list(i) = head(idxcell)
            head(idxcell) = i
        enddo

    end subroutine links

    subroutine forces_linkedlist(Ncell,mapa,head,list,x,f,w,L,rcut)
        use precision
        implicit none

        integer(dp), intent(in) :: Ncell,mapa(:),list(:),head(:)
        real(dp), intent(in)    :: x(:,:),L,rcut
        real(dp), intent(out)   :: f(:,:),w
        integer(dp)             :: i,j,idxcell,jcell,jcell0
        real(dp)                :: dx1,dx2,dx3,u,rij,nabor
        
        u = 0._dp
        w = 0._dp
        
        do idxcell = 1,Ncell
            i = head(idxcell)
            do while(i/=0_dp)
                j = list(i)
                do while(j/=0_dp)
                    call f_ij(L,j,i,rcut,x,f,w)
                    j = list(j)
                enddo
                jcell0 = 13_dp*(idxcell-1_dp)
                do nabor = 1,13
                    jcell = mapa(jcell0 + nabor)
                    j = head(jcell)
                    do while(j/=0)
                        call f_ij(L,j,i,rcut,x,f,w)
                        j = list(j)
                    enddo
                enddo
                i = list(i)
            enddo
        enddo
    end subroutine forces_linkedlist    

    subroutine velocity_verlet(r,v,f,dt,rcut,L,w)
        use precision
        implicit none

        integer(dp)                 :: i,natm
        real(dp), intent(inout)     :: r(:,:),v(:,:),f(:,:)
        real(dp), intent(in)        :: dt,rcut,L
        real(dp), intent(out)       :: w
        real(dp), allocatable       :: r1(:,:),v1(:,:),f1(:,:)

        natm = size(r(:,1))

        allocate(r1(natm,3),v1(natm,3),f1(natm,3))

        r1 = 0._dp;v1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            r1(i,:) = r(i,:) + v(i,:)*dt + f(i,:)*dt*dt*0.5_dp
        enddo
        call md_pbc(r1,L)
        call system_forces(r1,rcut,f1,L,w)
        do i = 1,size(r(:,1))
            v1(i,:) = v(i,:) + (f1(i,:)+f(i,:))*dt*0.5_dp
        enddo

        r = r1
        v = v1
        f = f1

    end subroutine velocity_verlet


    subroutine velocity_verlet_LL(mapa,head,list,r,v,f,dt,rcut,L,w,Mcell,Ncell)
        use precision
        implicit none

        integer(dp)                 :: i,natm
        real(dp), intent(inout)     :: r(:,:),v(:,:),f(:,:)
        integer(dp), intent(inout)  :: head(:),list(:)
        real(dp), intent(in)        :: dt,rcut,L
        integer(dp), intent(in)     :: Mcell,Ncell,mapa(:)
        real(dp), intent(out)       :: w
        real(dp), allocatable       :: r1(:,:),v1(:,:),f1(:,:)

        natm = size(r(:,1))

        allocate(r1(natm,3),v1(natm,3),f1(natm,3))

        r1 = 0._dp;v1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            r1(i,:) = r(i,:) + v(i,:)*dt + f(i,:)*dt*dt*0.5_dp
        enddo
        call md_pbc(r1,L)
        call links(Mcell,L,rcut,r1,head,list)
        call forces_linkedlist(Ncell,mapa,head,list,r1,f1,w,L,rcut)
        do i = 1,size(r(:,1))
            v1(i,:) = v(i,:) + (f1(i,:)+f(i,:))*dt*0.5_dp
        enddo

        r = r1
        v = v1
        f = f1

    end subroutine velocity_verlet_LL

    subroutine dU_LJ(r,dx1,dx2,dx3,rcut,index,L,delta_U)
        use precision
        implicit none

        real(dp), intent(in)    :: r(:,:),dx1,dx2,dx3,rcut,L
        integer(dp), intent(in) :: index
        real(dp), intent(out)   :: delta_U
        integer(dp)             :: i,j,natm
        real(dp)                :: Uscale,U_old,U_new,Ui_old,Ui_new
        real(dp), allocatable   :: r_new(:,:)

        !call U_LJ(rcut,Uscale)
        Uscale = 0._dp

        U_new = 0._dp; U_old = 0._dp; Ui_old = 0._dp; Ui_new = 0._dp

        natm = int(size(r(:,1)))
        allocate(r_new(natm,3))
        r_new = 0._dp

        r_new = r

        r_new(index,1) = r_new(index,1) + dx1
        r_new(index,2) = r_new(index,2) + dx2
        r_new(index,3) = r_new(index,3) + dx3

        call md_pbc(r_new,L)

        do i = 1,natm
            if (i /= index) then
                call U_ij(L,index,i,rcut,r,Uscale,Ui_old)
                call U_ij(L,index,i,rcut,r_new,Uscale,Ui_new)
                U_old = U_old + Ui_old; U_new = U_new + Ui_new
            endif
        enddo
        
        delta_U = (U_new - U_old)

    end subroutine dU_LJ    

    subroutine mc_step_LJ(r,dx,T,dU_total,seed,L,rcut)
        use precision,prs=>sp
        use mtmod,only: sgrnd,grnd
        implicit none

        real(dp), intent(inout)     :: r(:,:)
        real(dp), intent(in)        :: T,dx,L,rcut
        real(dp), intent(out)       :: dU_total
        integer(4),intent(inout)    :: seed
        integer(dp)                 :: index,i,j,k,nn
        real(dp)                    :: delta_U,Bolt
        real(dp)                    :: dx1,dx2,dx3,ran_u
        

        call sgrnd(seed)

        dU_total = 0._dp

        nn = size(r(:,1))

        do k = 1,nn
            delta_U = 0._dp
            index = floor(grnd()*real(nn,dp))+1_dp
            
            dx1 = dx*real(grnd()-0.5_dp)
            dx2 = dx*real(grnd()-0.5_dp)
            dx3 = dx*real(grnd()-0.5_dp)

            call dU_LJ(r,dx1,dx2,dx3,rcut,index,L,delta_U)
            if (delta_U <= 0._dp) then
                r(index,1) = r(index,1) + dx1
                r(index,2) = r(index,2) + dx2
                r(index,3) = r(index,3) + dx3
                call md_pbc(r,L)
                dU_total = dU_total + delta_U
            else if (delta_U > 0._dp) then
                ran_u = grnd()
                if (T == 0._dp) then
                    Bolt = 0._dp
                else
                    Bolt = exp(-delta_U/T)
                endif
                if (Bolt >= ran_u) then
                    r(index,1) = r(index,1) + dx1
                    r(index,2) = r(index,2) + dx2
                    r(index,3) = r(index,3) + dx3
                    call md_pbc(r,L)
                    dU_total = dU_total + delta_U
                endif
            endif
        enddo
        seed = abs(int(grnd()*9000))+150_dp
        dU_total = dU_total/real(nn,dp)
        
    end subroutine mc_step_LJ

    subroutine mc_step_LJ_nopbc(r,r0,dx,T,dU_total,seed,L,rcut)
        use precision,prs=>sp
        use mtmod,only: sgrnd,grnd
        implicit none

        real(dp), intent(inout)     :: r(:,:),r0(:,:)
        real(dp), intent(in)        :: T,dx,L,rcut
        real(dp), intent(out)       :: dU_total
        integer(4),intent(inout)    :: seed
        integer(dp)                 :: index,i,j,k,nn
        real(dp)                    :: delta_U,Bolt
        real(dp)                    :: dx1,dx2,dx3,ran_u
        

        call sgrnd(seed)

        dU_total = 0._dp

        nn = size(r(:,1))

        do k = 1,nn
            delta_U = 0._dp
            index = floor(grnd()*real(nn,dp))+1_dp
            
            dx1 = dx*real(grnd()-0.5_dp)
            dx2 = dx*real(grnd()-0.5_dp)
            dx3 = dx*real(grnd()-0.5_dp)

            call dU_LJ(r,dx1,dx2,dx3,rcut,index,L,delta_U)
            if (delta_U <= 0._dp) then
                r(index,1) = r(index,1) + dx1
                r(index,2) = r(index,2) + dx2
                r(index,3) = r(index,3) + dx3
                r0(index,1) = r0(index,1) + dx1
                r0(index,2) = r0(index,2) + dx2
                r0(index,3) = r0(index,3) + dx3
                call md_pbc(r,L)
                dU_total = dU_total + delta_U
            else if (delta_U > 0._dp) then
                ran_u = grnd()
                if (T == 0._dp) then
                    Bolt = 0._dp
                else
                    Bolt = exp(-delta_U/T)
                endif
                if (Bolt >= ran_u) then
                    r(index,1) = r(index,1) + dx1
                    r(index,2) = r(index,2) + dx2
                    r(index,3) = r(index,3) + dx3
                    r0(index,1) = r0(index,1) + dx1
                    r0(index,2) = r0(index,2) + dx2
                    r0(index,3) = r0(index,3) + dx3
                    call md_pbc(r,L)
                    dU_total = dU_total + delta_U
                endif
            endif
        enddo
        seed = abs(int(grnd()*9000))+150_dp
        dU_total = dU_total/real(nn,dp)
        
    end subroutine mc_step_LJ_nopbc

    subroutine mc_step_LJ_test(r,dx,T,seed,L,rcut,ratio)
        use precision,prs=>sp
        use mtmod,only: sgrnd,grnd
        implicit none

        real(dp), intent(inout)     :: r(:,:)
        real(dp), intent(in)        :: T,dx,L,rcut
        integer(sp), intent(inout)  :: seed
        real(dp), intent(out)       :: ratio
        integer(dp)                 :: index,i,j,k,nn
        real(dp)                    :: delta_U,Bolt,U
        real(dp)                    :: dx1,dx2,dx3,ran_u                

        call sgrnd(seed)

        nn = size(r(:,1))
        U = 0._dp
        ratio = 0._dp

        do k = 1,nn
            delta_U = 0._dp
            index = floor(grnd()*real(nn,dp))+1_dp
            
            dx1 = dx*real(grnd()-0.5_dp)
            dx2 = dx*real(grnd()-0.5_dp)
            dx3 = dx*real(grnd()-0.5_dp)

            call dU_LJ(r,dx1,dx2,dx3,rcut,index,L,delta_U)
            if (delta_U <= 0._dp) then
                ratio = ratio + 1._dp
                r(index,1) = r(index,1) + dx1
                r(index,2) = r(index,2) + dx2
                r(index,3) = r(index,3) + dx3

                call md_pbc(r,L)
                U = U + delta_U
            else
                ran_u = grnd()
                if (T == 0._dp) then
                    Bolt = 0._dp
                else
                    Bolt = exp(-delta_U/T)
                endif
                if (Bolt >= ran_u) then
                    ratio = ratio + 1._dp
                    r(index,1) = r(index,1) + dx1
                    r(index,2) = r(index,2) + dx2
                    r(index,3) = r(index,3) + dx3
                    call md_pbc(r,L)
                    U = U + delta_U
                endif
            endif
        enddo
        seed = abs(int(grnd()*9000))+150_dp

        ratio = ratio/real(nn,dp)

    end subroutine mc_step_LJ_test


    subroutine set_dx(r,Tref,seed,L,rcut,dx)
        use precision
        implicit none

        real(dp), intent(in)    :: r(:,:),Tref,L,rcut
        real(dp), allocatable   :: r_n(:,:)
        integer(sp), intent(in) :: seed
        real(dp), intent(out)   :: dx
        integer(sp)             :: seed_n
        real(dp)                :: ratio

        allocate(r_n(size(r),3))
        r_n = r
        dx = 1._dp

        seed_n = seed
        do while(ratio > 0.55_dp .or. ratio < 0.45_dp)
            call mc_step_LJ_test(r_n,dx,Tref,seed_n,L,rcut,ratio)
            if (ratio > 0.55_dp) then
                dx = dx*1.05_dp
            else if (ratio < 0.45_dp) then
                dx = dx*0.95_dp
            endif
        enddo

    end subroutine set_dx

    subroutine system_U(r,rcut,U,L)
        use precision
        implicit none

        real(dp), intent(in)    :: r(:,:)
        real(dp), intent(in)    :: rcut,L
        real(dp), intent(out)   :: U
        integer(dp)             :: i,j
        real(dp)                :: rij,dx1,dx2,dx3,Uscale,dU

        U = 0._dp; Uscale = 0._dp

        !call U_LJ(rcut,Uscale)

        do j = 2,size(r(:,1))
            do i = 1,j-1
                dU = 0._dp
                call U_ij(L,j,i,rcut,r,Uscale,dU)
                U = U + dU
            enddo
        enddo

        U = U/real(size(r(:,1)),dp)

    end subroutine system_U

    subroutine potencial_linkedlist(Ncell,mapa,head,list,x,U,L,rcut)
        use precision
        implicit none

        integer(dp), intent(in) :: Ncell,mapa(:),list(:),head(:)
        real(dp), intent(in)    :: x(:,:),L,rcut
        real(dp), intent(out)   :: U
        integer(dp)             :: i,j,idxcell,jcell,jcell0,nabor
        real(dp)                :: dx,dy,dz,rij,Uscale,dU

        Uscale = 0._dp; U = 0._dp

        !call U_LJ(rcut,Uscale)

        do idxcell = 1,Ncell
            i = head(idxcell)
            do while(i/=0_dp)
                j = list(i)
                do while(j/=0_dp)
                    call U_ij(L,j,i,rcut,x,Uscale,dU)
                    U = U + dU
                    j = list(j)
                enddo
                jcell0 = 13_dp*(idxcell-1_dp)
                do nabor = 1,13
                    jcell = mapa(jcell0 + nabor)
                    j = head(jcell)
                    do while(j/=0)
                        call U_ij(L,j,i,rcut,x,Uscale,dU)
                        U = U + dU
                        j = list(j)
                    enddo
                enddo
                i = list(i)
            enddo
        enddo

        U = U/real(size(x(:,1)),dp)
    end subroutine potencial_linkedlist


    subroutine system_Ek(v,Ek)
        use precision
        implicit none

        real(dp), intent(in)    :: v(:,:)
        real(dp), intent(out)   :: Ek
        integer(dp)             :: i

        Ek = 0._dp

        do i = 1,size(v(:,1))
            Ek = Ek +  v(i,1)*v(i,1) + v(i,2)*v(i,2) + v(i,3)*v(i,3)
        enddo

        Ek = Ek*0.5_dp*(1._dp/real(size(v(:,1)),dp))

    end subroutine system_Ek

    subroutine v_rescaling(Tref,v,N)
        use precision
        implicit none
        real(dp), intent(in)    :: Tref
        integer(dp), intent(in) :: N
        real(dp)                :: Temp
        real(dp), intent(inout) :: v(N,3)
        real(dp)                :: vi2
        integer(dp)             :: i

        vi2 = 0._dp

        do i = 1,N
            vi2 = vi2 + (v(i,1)*v(i,1)+v(i,2)*v(i,2)+v(i,3)*v(i,3))
        enddo

        Temp = vi2/(3._dp*real(N,dp))

        do i = 1,N
            v(i,1) = v(i,1)*sqrt(Tref/Temp)
            v(i,2) = v(i,2)*sqrt(Tref/Temp)
            v(i,3) = v(i,3)*sqrt(Tref/Temp)
        enddo

    end subroutine v_rescaling

    subroutine Temperature(v,N,rho,KbT)
        use precision
        implicit none

        integer(dp), intent(in) :: N
        integer(dp)             :: i
        real(dp), intent(in)    :: v(N,3)
        real(dp), intent(in)    :: rho
        real(dp), intent(out)   :: Kbt
        real(dp)                :: v2

        v2 = 0._dp

        do i = 1,N
            v2 = v2 + (v(i,1)*v(i,1)+v(i,2)*v(i,2)+v(i,3)*v(i,3))
        enddo

        KbT = v2/(3._dp*real(N,dp))

    end subroutine Temperature

    subroutine mc_speed(v,N,cm_v)
        use precision
        implicit none

        integer(dp), intent(in) :: N
        real(dp), intent(in)    :: v(N,3)
        real(dp), intent(out)   :: cm_v
        real(dp)                :: v_cm(3)            


        v_cm(1) = sum(v(:,1))/real(N,dp)
        v_cm(2) = sum(v(:,2))/real(N,dp)
        v_cm(3) = sum(v(:,3))/real(N,dp)

        cm_v = sqrt(v_cm(1)*v_cm(1)+v_cm(2)*v_cm(2)+v_cm(3)*v_cm(3))

    end subroutine mc_speed

    subroutine histogram_3d(x,n_bins,counts,bins)
        use precision
        implicit none

        integer(dp)                             :: i,j,k,N
        integer(dp), intent(in)                 :: n_bins
        real(dp), intent(in)                    :: x(:,:)
        real(dp), allocatable, intent(out)      :: bins(:,:)
        integer(dp), allocatable, intent(out)   :: counts(:,:)
        real(dp)                                :: min1,min2,min3,max1,max2,max3
        real(dp), allocatable                   :: x_abs(:,:)

        N = size(x(:,1))

        allocate(x_abs(N,3))

        x_abs = abs(x)

        allocate(bins(n_bins,3),counts(n_bins,3))

        bins = 0._dp;counts = 0._dp

        min1 = minval(x_abs(:,1));min2 = minval(x_abs(:,2));min3 = minval(x_abs(:,3))
        max1 = maxval(x_abs(:,1));max2 = maxval(x_abs(:,2));max3 = maxval(x_abs(:,3))

        do i = 1,n_bins
            bins(i,1) = min1 + (max1 - min1)*real(i-1,dp)/real(n_bins,dp) + (max1 - min1)*0.5_dp*(1._dp/real(n_bins,dp))
            bins(i,2) = min2 + (max2 - min2)*real(i-1,dp)/real(n_bins,dp) + (max2 - min2)*0.5_dp*(1._dp/real(n_bins,dp))
            bins(i,3) = min3 + (max3 - min3)*real(i-1,dp)/real(n_bins,dp) + (max3 - min3)*0.5_dp*(1._dp/real(n_bins,dp))
        enddo

        do i = 1,N
            do j = 1,n_bins
            if (j /= n_bins) then
                if (x_abs(i,1) >= bins(j,1) .and. x_abs(i,1) < bins(j+1,1)) then
                    counts(j,1) = counts(j,1) + 1_dp
                endif
                if (x_abs(i,2) >= bins(j,2) .and. x_abs(i,2) < bins(j+1,2)) then
                    counts(j,2) = counts(j,2) + 1_dp
                endif
                if (x_abs(i,3) >= bins(j,3) .and. x_abs(i,3) < bins(j+1,3)) then
                    counts(j,3) = counts(j,3) + 1_dp
                endif
            else if (j == n_bins) then
                if (x_abs(i,1) >= bins(j-1,1) .and. x_abs(i,1) < bins(j,1)) then
                    counts(j,1) = counts(j,1) + 1_dp
                endif
                if (x_abs(i,2) >= bins(j-1,2) .and. x_abs(i,2) < bins(j,2)) then
                    counts(j,2) = counts(j,2) + 1_dp
                endif
                if (x_abs(i,3) >= bins(j-1,3) .and. x_abs(i,3) < bins(j,3)) then
                    counts(j,3) = counts(j,3) + 1_dp
                endif
            endif
        enddo
    enddo

    end subroutine histogram_3d

    subroutine gr(N,r,dr,g,nbins,L,rho)
        use precision
        implicit none

        integer(dp), intent(in) :: N,nbins
        real(dp), intent(in)    :: L,rho
        real(dp), intent(in)    :: r(:,:)
        real(dp), intent(out)   :: dr(:),g(:)
        real(dp), parameter     :: pi=4._dp*atan(1._dp)
        integer(dp)             :: i,j,k
        real(dp)                :: r_max,rij2,rij
        real(dp)                :: dx1,dx2,dx3,bin

        r_max = 0._dp

        dr = 0._dp; g = 0._dp

        do j = 2,N
            do i = 1,j-1
                dx1 = r(j,1)-r(i,1)
                dx2 = r(j,2)-r(i,2)
                dx3 = r(j,3)-r(i,3)
                dx1 = dx1 - L*anint(dx1/L)
                dx2 = dx2 - L*anint(dx2/L)
                dx3 = dx3 - L*anint(dx3/L)
                rij2 = dx1*dx1 + dx2*dx2 + dx3*dx3
                if (sqrt(rij2) > r_max) then
                    r_max = sqrt(rij2)
                endif
            enddo
        enddo

        bin = r_max/(real(nbins,dp))

        do i = 1,nbins-1
            dr(i) = r_max*(real(i-1,dp)/real(nbins,dp))
        enddo

        do i = 1,N
            do j = 1,N
                if (i /= j) then
                    dx1 = r(j,1)-r(i,1)
                    dx2 = r(j,2)-r(i,2)
                    dx3 = r(j,3)-r(i,3)
                    dx1 = dx1 - L*anint(dx1/L)
                    dx2 = dx2 - L*anint(dx2/L)
                    dx3 = dx3 - L*anint(dx3/L)
                    rij = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3)
                    do k = 1,nbins
                        if (k < nbins) then
                            if (rij >= dr(k) .and. rij < dr(k+1)) then
                                g(k) = g(k) + 1._dp
                            endif
                        elseif (k == nbins) then
                            if (rij >= dr(k-1) .and. rij <= dr(k)) then
                                g(k) = g(k) + 1._dp
                            endif
                        endif
                    enddo
                endif
                enddo
            enddo
            
            do i = 1,nbins
                g(i) = (g(i)*3._dp)/(N*4._dp*rho*pi*((dr(i)+bin)**3 - dr(i)**3))
            enddo

    end subroutine gr

    subroutine parametro_cristalino(r,natm,L,S)
        use precision
        implicit none

        real(dp), intent(in)    :: r(:,:)
        integer(dp), intent(in) :: natm
        real(dp), intent(in)    :: L
        real(dp), intent(out)   :: S
        real(dp)                :: a
        integer(dp)             :: i
        real(dp), allocatable   :: k(:)
        real(dp), parameter     :: pi=4._dp*atan(1._dp)
        real(dp)                :: re,im,mod2

        re = 0._dp;im = 0._dp
        allocate(k(3))

        a = L/(real(natm,dp)*0.25_dp)**(1._dp/3._dp)

        k = (2._dp*pi/a)*[-1._dp,1._dp,-1._dp]

        do i = 1,natm
            re = re + cos(dot_product(k,r(i,:)))
            im = im + sin(dot_product(k,r(i,:)))
        enddo
        
        mod2 = re*re + im*im
        S = mod2/(real(natm*natm,dp))


    end subroutine parametro_cristalino

    subroutine MSD(natm,nmsd,r,Wxx,Wyy,Wzz,sumWxx,sumWyy,sumWzz,dindex,enabled)
        use precision
        implicit none

        integer(dp), save           :: idif
        integer(dp)                 :: last,i,j,inxt
        integer(dp), intent(in)     :: nmsd,natm,enabled
        real(dp), intent(in)        :: r(:,:)
        real(dp), intent(inout)     :: Wxx(:,:),Wyy(:,:),Wzz(:,:)
        real(dp), intent(inout)     :: sumWxx(:),sumWyy(:),sumWzz(:)
        integer(dp), intent(inout)  :: dindex(:)

        if (enabled == 1_dp) idif = 0_dp

        idif = idif + 1_dp
        last = mod(idif-1,nmsd) + 1_dp

        do i = 1,natm
            Wxx(i,last) = r(i,1)
            Wyy(i,last) = r(i,2)
            Wzz(i,last) = r(i,3)
        enddo

        if (mod(idif,10_dp) == 0_dp) then
            if (idif > nmsd) then
                do j = 1,nmsd
                    inxt = mod(idif-j,nmsd) + 1_dp
                    do i = 1,natm
                        sumWxx(j) = sumWxx(j) + (Wxx(i,last)-Wxx(i,inxt))**2
                        sumWyy(j) = sumWyy(j) + (Wyy(i,last)-Wyy(i,inxt))**2
                        sumWzz(j) = sumWzz(j) + (Wzz(i,last)-Wzz(i,inxt))**2
                    enddo
                    dindex(j) = dindex(j) + 1_dp
                enddo
            endif
        endif
        
    end subroutine MSD

    subroutine velocity_verlet_nopbc(r,v,f,dt,rcut,L,w,r0)
        use precision
        implicit none

        integer(dp)             :: i,natm
        real(dp), intent(inout) :: r(:,:),v(:,:),f(:,:)
        real(dp), intent(in)    :: dt,rcut,L
        real(dp), intent(inout) :: r0(:,:)
        real(dp), intent(out)   :: w
        real(dp), allocatable   :: r1(:,:),v1(:,:),f1(:,:)

        natm = size(r(:,1))

        allocate(r1(natm,3),v1(natm,3),f1(natm,3))

        r1 = 0._dp;v1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            r1(i,:) = r(i,:) + v(i,:)*dt + f(i,:)*dt*dt*0.5_dp
            r0(i,:) = r0(i,:) + v(i,:)*dt + f(i,:)*dt*dt*0.5_dp
        enddo
        call md_pbc(r1,L)
        call system_forces(r1,rcut,f1,L,w)
        do i = 1,size(r(:,1))
            v1(i,:) = v(i,:) + (f1(i,:)+f(i,:))*dt*0.5_dp
        enddo

        r = r1
        v = v1
        f = f1

    end subroutine velocity_verlet_nopbc

    subroutine velocity_verlet_LL_nopbc(mapa,head,list,r,v,f,dt,rcut,L,w,Mcell,Ncell,r0)
        use precision
        implicit none

        integer(dp)                         :: i,natm
        real(dp), intent(inout)             :: r(:,:),v(:,:),f(:,:)
        real(dp), intent(out), allocatable  :: r0(:,:)
        integer(dp), intent(inout)          :: head(:),list(:)
        real(dp), intent(in)                :: dt,rcut,L
        integer(dp), intent(in)             :: Mcell,Ncell,mapa(:)
        real(dp), intent(out)               :: w
        real(dp), allocatable               :: r1(:,:),v1(:,:),f1(:,:)

        natm = size(r(:,1))

        allocate(r1(natm,3),v1(natm,3),f1(natm,3),r0(natm,3))

        r1 = 0._dp;v1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            r1(i,:) = r(i,:) + v(i,:)*dt + f(i,:)*dt*dt*0.5_dp
            r0(i,:) = r0(i,:) + v(i,:)*dt + f(i,:)*dt*dt*0.5_dp
        enddo
        call md_pbc(r1,L)
        call links(Mcell,L,rcut,r1,head,list)
        call forces_linkedlist(Ncell,mapa,head,list,r1,f1,w,L,rcut)
        do i = 1,size(r(:,1))
            v1(i,:) = v(i,:) + (f1(i,:)+f(i,:))*dt*0.5_dp
        enddo

        r = r1
        v = v1
        f = f1

    end subroutine velocity_verlet_LL_nopbc

    subroutine Xi(D0,dt,idum,x)
        use precision
        use rangens, only: gaussdev
        implicit none

        integer(sp), intent(in) :: idum
        real(dp), intent(in)    :: D0,dt
        real(dp), intent(out)   :: x
        real(dp)                :: xn

        x = gaussdev(idum)
        xn = x
        x = sqrt(2._dp*D0*dt)*x
    end subroutine

    subroutine brownian_r(r,f,dt,rcut,L,w,KbT,D0,idum)
        use precision
        implicit none

        real(dp), parameter         :: pi=4._dp*atan(1._dp)
        integer(sp), intent(in)     :: idum
        integer(dp)                 :: i,natm
        real(dp), intent(inout)     :: r(:,:),f(:,:)
        real(dp), intent(in)        :: dt,rcut,L,D0,KbT
        real(dp), intent(out)       :: w
        real(dp), allocatable       :: r1(:,:),f1(:,:)
        real(dp)                    :: x_i

        natm = size(r(:,1))

        allocate(r1(natm,3),f1(natm,3))

        r1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            call Xi(D0,dt,idum,x_i)
            r1(i,:) = r(i,:) + f(i,:)*D0*dt/KbT + x_i
        enddo
        call md_pbc(r1,L)
        call system_forces(r1,rcut,f1,L,w)

        r = r1
        f = f1
    end subroutine brownian_r

    subroutine brownian_r_LL(mapa,head,list,Mcell,Ncell,r,f,dt,rcut,L,w,KbT,D0,idum)
        use precision
        use mtmod ,only: sgrnd,grnd
        implicit none

        real(dp), parameter         :: pi=4._dp*atan(1._dp)
        integer(sp), intent(inout)  :: idum
        integer(dp)                 :: i,natm
        real(dp), intent(inout)     :: r(:,:),f(:,:)
        integer(dp), intent(inout)  :: head(:),list(:)
        integer(dp), intent(in)     :: Mcell,Ncell,mapa(:)
        real(dp), intent(in)        :: dt,rcut,L,D0,KbT
        real(dp), intent(out)       :: w
        real(dp), allocatable       :: r1(:,:),f1(:,:)
        real(dp)                    :: x_i

        natm = size(r(:,1))
        w = 0._dp

        allocate(r1(natm,3),f1(natm,3))

        r1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            call Xi(D0,dt,idum,x_i)
            r1(i,:) = r(i,:) + f(i,:)*D0*dt/KbT + x_i
        enddo
        call md_pbc(r1,L)
        call links(Mcell,L,rcut,r1,head,list)
        call forces_linkedlist(Ncell,mapa,head,list,r1,f1,w,L,rcut)
        r = r1
        f = f1

        call sgrnd(idum)

        idum = abs(int(grnd()*9000))+150
    end subroutine brownian_r_LL

    subroutine brownian_r_nopbc(r,r_nopbc,f,dt,rcut,L,w,KbT,D0,idum)
        use precision
        implicit none

        real(dp), parameter         :: pi=4._dp*atan(1._dp)
        integer(sp), intent(in)     :: idum
        integer(dp)                 :: i,natm
        real(dp), intent(inout)     :: r(:,:),f(:,:),r_nopbc(:,:)
        real(dp), intent(in)        :: dt,rcut,L,D0,KbT
        real(dp), intent(out)       :: w
        real(dp), allocatable       :: r1(:,:),f1(:,:)
        real(dp)                    :: x_i

        natm = size(r(:,1))

        allocate(r1(natm,3),f1(natm,3))

        r1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            call Xi(D0,dt,idum,x_i)
            r1(i,:) = r(i,:) + f(i,:)*D0*dt/KbT + x_i
            r_nopbc(i,:) = r_nopbc(i,:) + f(i,:)*D0*dt/KbT + x_i
        enddo
        call md_pbc(r1,L)
        call system_forces(r1,rcut,f1,L,w)

        r = r1
        f = f1
    end subroutine brownian_r_nopbc

   
    subroutine brownian_r_LL_nopbc(mapa,head,list,Mcell,Ncell,r,r_nopbc,f,dt,rcut,L,w,KbT,D0,idum)
        use precision
        use mtmod ,only: sgrnd,grnd
        implicit none

        real(dp), parameter         :: pi=4._dp*atan(1._dp)
        integer(sp), intent(inout)     :: idum
        integer(dp)                 :: i,natm
        real(dp), intent(inout)     :: r(:,:),f(:,:),r_nopbc(:,:)
        integer(dp), intent(inout)  :: head(:),list(:)
        integer(dp), intent(in)     :: Mcell,Ncell,mapa(:)
        real(dp), intent(in)        :: dt,rcut,L,D0,KbT
        real(dp), intent(out)       :: w
        real(dp), allocatable       :: r1(:,:),f1(:,:)
        real(dp)                    :: x1_i,x2_i,x3_i

        natm = size(r(:,1))

        allocate(r1(natm,3),f1(natm,3))

        r1 = 0._dp;f1 = 0._dp

        do i = 1,size(r(:,1))
            call Xi(D0,dt,idum,x1_i)
            call Xi(D0,dt,idum,x2_i)
            call Xi(D0,dt,idum,x3_i)
            r1(i,1) = r(i,1) + f(i,1)*D0*dt/KbT + x1_i
            r1(i,2) = r(i,2) + f(i,2)*D0*dt/KbT + x2_i
            r1(i,3) = r(i,3) + f(i,3)*D0*dt/KbT + x3_i
            r_nopbc(i,1) = r_nopbc(i,1) + f(i,1)*D0*dt/KbT + x1_i
            r_nopbc(i,2) = r_nopbc(i,2) + f(i,2)*D0*dt/KbT + x2_i
            r_nopbc(i,3) = r_nopbc(i,3) + f(i,3)*D0*dt/KbT + x3_i
        enddo
        call md_pbc(r1,L)
        call links(Mcell,L,rcut,r1,head,list)
        call forces_linkedlist(Ncell,mapa,head,list,r1,f1,w,L,rcut)

        r = r1
        f = f1

        call sgrnd(idum)

        idum = abs(int(grnd()*9000))+150
    end subroutine brownian_r_LL_nopbc
    
    subroutine binomial(p,n,seed,x)
        use precision
        use ran2mod, only: ran2
        implicit none

        real(dp), intent(in)        :: p
        integer(dp), intent(in)     :: n
        integer(sp), intent(in)     :: seed
        integer(dp), intent(out)    :: x
        integer(dp)                 :: i
        real(dp)                    :: u

        x = 0._dp

        do i = 1,n
            u = ran2(seed)
            if (u <= p) then
                x = x + 1_dp
            endif
        enddo

    end subroutine binomial


    subroutine histograma(x,nbins,bins,counts)
        use precision
        implicit none

        integer(dp)                         :: i,j,n_dat
        integer(dp), intent(in)             :: nbins
        real(dp), intent(in)                :: x(:)
        real(dp), intent(out), allocatable  :: bins(:),counts(:)
        real(dp)                            :: binwidth,min,max

        allocate(bins(nbins),counts(nbins))
        bins = 0._dp; counts = 0._dp

        n_dat = size(x)
        
        min = minval(x)
        max = maxval(x)

        binwidth = (max - min)/real(n_dat,dp)

        do i = 1,nbins
            bins(i) = min + (max - min)*real(i,dp)/real(nbins,dp)
        enddo

        do i = 1,nbins-1
            do j = 1,n_dat
                if (bins(i) <= x(j) .and. bins(i+1) > x(j)) then
                    counts(i) = counts(i) + 1._dp
                endif
            enddo
        enddo

        bins = bins + binwidth*0.5_dp
        counts = counts/(binwidth*real(n_dat,dp))


    end subroutine histograma

    function p_gaussian(x,mu,sigma)
        use precision
        implicit none

        real(dp), intent(in)    :: mu,sigma,x
        real(dp)                :: p_gaussian
        real(dp), parameter     :: pi=4._dp*atan(1._dp)
        

        p_gaussian = exp((-(x-mu)**2)/(2._dp*sigma*sigma))/sqrt(2._dp*pi*sigma*sigma)

    end function


    function OU_P(x,lamda)
        use precision
        use rangens, only: gaussdev
        implicit none

        real(dp)            :: x,lamda
        real(dp)            :: OU_P
        real(dp)            :: gw
        integer(sp), save   :: seed

        seed = 34521_sp

        OU_P = -lamda*x + lamda*gaussdev(seed)
        
    end function

    subroutine linear_regression(x,y,a,b)
        use precision
        implicit none

        real(dp), intent(in)    :: x(:),y(:)
        real(dp), intent(out)   :: a,b
        integer(dp)             :: i,N
        real(dp)                :: nn,sum_y,sum_x,sum_xy,sum_x2

        N = int(size(x),dp)
        nn = real(N,dp)

        sum_x = sum(x);sum_y = sum(y)
        sum_xy = 0._dp;sum_x2 = 0._dp

        do i = 1,N
            sum_xy = sum_xy + (x(i)*y(i))
            sum_x2 = sum_x2 + x(i)**2
        enddo
    
        b = ((sum_y*sum_x2)-(sum_x*sum_xy))/(nn*sum_x2-sum_x*sum_x)
        a = ((nn*sum_xy)-(sum_x*sum_y))/(nn*sum_x2-sum_x*sum_x)

    end subroutine linear_regression

end module miscelaneous