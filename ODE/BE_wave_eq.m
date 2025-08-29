function [U,V]=BE_wave_eq(Dati,x,dependence,Matrix,list_U ,list_V,front,t_hat)
    
    T = 1; %T is equal 1 indeed RK is applied in K_i_hat so the maximum time is 1
    dt = Dati.dt; 
    % initialization
    t = 0:dt:T; % Vettore temporale
    num_steps = length(t);
    u0=zeros(length(x),1);
    v0=u0;
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    % initial condition
    if sum((t_hat==0))>=2
        [u0,v0] = impose_ic(Dati,front,x,Matrix);
    else 
        [u0,v0] = impose_continuity1(dependence,list_U,list_V,x,Dati);
    end

    Y = [u0; v0]; 
    %Pe = Dati.rho*mean(u0)*Dati.h/(2*Dati.mu)
    % preallocation of the buffer
    U = zeros(length(u0), num_steps);
    V = zeros(length(v0), num_steps);
    dimx = length(u0);
    % store the initial condition
    U(:,1) = u0;
    V(:,1) = v0;

    for n = 1:num_steps-1
        % actual time
        ti = t(n);
        tii = t(n+1);
        Yi = Y(:,n);
        
        
        Matrices = compute_matrix1(Dati,Matrix,front,x,tii);
        M = Matrices.M;
        D = Matrices.D;
        W = Matrices.W;
        eps =  dt*0.5*Dati.h  * sqrt((D*Yi)'*W*(D*Yi));
        eps_min = 1e-4;
        eps_max = 1e-2;
        eps = max(min(eps, eps_max), eps_min);

        %eps = 0;
        Sigma = Matrices.Sigma;
        I = numerical_flux(front,Dati,x);
        A = M + dt*(Sigma-I);
        
        Y0 = zeros(length(Yi),1);
        Y1 = zeros(length(Yi),1);
        if sum(x==x0)>0
            Ybc = impose_bc(x,tii,Dati,front);
            Y0(1) = Ybc(1);
            Y0(length(u0)) = Yi(dimx);
            Y0(length(u0)+1) = Ybc(dimx+1);
            Y0(end) = Yi(end);
        elseif sum(x==x1)>0
            Ybc = impose_bc(x,tii,Dati,front);
            Y0(1) = Yi(1);
            Y0(length(u0)) = Ybc(dimx);
            Y0(length(u0)+1) = Yi(dimx+1);
            Y0(end) = Ybc(end);
        else
            Y0(1) = Yi(1);
            Y0(length(u0)) = Yi(dimx);
            Y0(length(u0)+1) = Yi(dimx+1);
            Y0(end) = Yi(end);
        end

        idx_fissi = [1, dimx, dimx+1, 2*dimx];
        idx_tutti = 1:2*dimx;
        idx_liberi = setdiff(idx_tutti, idx_fissi);
        A_red = A(idx_liberi,idx_liberi);
        
        if ~isempty(idx_liberi)
            b = M*Yi - A*Y0;
            b_red = b(idx_liberi);
            Y1(idx_liberi) = A_red\b_red;
        end
        % update the solution
        Y(:,n+1) = Y1 + Y0;
        


        % split u and v
        U(:,n+1) = Y(1:length(u0),n+1);
        V(:,n+1) = Y(length(u0)+1:end, n+1);

    end


end