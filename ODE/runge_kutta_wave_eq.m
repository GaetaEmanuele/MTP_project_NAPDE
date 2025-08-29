function [U,V]=runge_kutta_wave_eq(Dati,x,dependence,Matrix,list_U ,list_V,front,t_hat)
    
    T = 1; %T is equal 1 indeed RK is applied in K_i_hat so the maximum time is 1
    dt = Dati.dt; 
    % initialization
    t = 0:Dati.dt:T; % Vettore temporale
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
    
    % store the initial condition
    U(:,1) = u0;
    V(:,1) = v0;

    for n = 1:num_steps-1
        % actual time
        ti = t(n);
        Yi = Y(:,n);
        
        
        % evalution of the RK43 coffiecent
        k1 = rk4_step(Yi, ti,x,Dati,Matrix,front);

        k2 = rk4_step(Yi + dt * k1/2, ti + (1/2)* dt,x,Dati,Matrix,front);

        k3 = rk4_step(Yi +(1/2) * dt * k2 , ti + (1/2) * dt,x,Dati,Matrix,front);

        k4 = rk4_step(Yi +dt*k3, ti + dt,x,Dati,Matrix,front);

        
        % update the solution
        Y(:,n+1) = Yi + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
        
        Y(1,n+1) = Y(1,n);
        Y(length(u0),n+1) = Y(length(u0),n);
        Y(length(u0)+1,n+1) = Y(length(u0)+1,n);
        Y(end,n+1) = Y(end,n);

        if sum(abs(x-x0)<1e-3) || sum(abs(x-x1)<1e-3)
            Y_ = impose_bc1(x,t(n+1),Dati,front);
            idx = find(Y_ ~= 0);
            Y(idx,n+1) = Y_(idx);
        end

        % split u and v
        U(:,n+1) = Y(1:length(u0),n+1);
        V(:,n+1) = Y(length(u0)+1:end, n+1);
    end
   
end
