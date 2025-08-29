function [U,V]=LDDRK_wave_eq(Dati,x,dependence,Matrix,list_U ,list_V,front,t_hat)
    
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

    a = [0.0, -0.737101392796, -1.634740794341, -0.744739003780, -1.469897351522, -2.813971388035];
    b = [0.032918605146, 0.823256998200, 0.381530948900, 0.200092213184, 1.718581042715, 0.270000000000];
    c = [0.0, 0.032918605146, 0.249351723343, 0.466911705055, 0.582030414044, 0.847252983783];


    for n = 1:num_steps-1
        % actual time
        ti = t(n);
        Yi = Y(:,n);
        w = zeros(size(Yi));
        Yj = Yi;
        
        for j=1:6
            dw_i = rk4_step(Yj,ti + c(j)*dt,x,Dati,Matrix,front);
            w = a(j)*w + dt*dw_i;
            Yj = Yj + b(j)*w;
        end
        
        Y(:,n+1) = Yj;


        Y(1,n+1) = Y(1,n);
        Y(length(u0),n+1) = Y(length(u0),n);
        Y(length(u0)+1,n+1) = Y(length(u0)+1,n);
        Y(end,n+1) = Y(end,n);
        
        if sum(abs(x-x0)<1e-3) || sum(abs(x-x1)<1e-3)
            Y_ = impose_bc(x,t(n+1),Dati,front);
            idx = find(Y_ ~= 0);
            Y(idx,n+1) = Y_(idx);

        end
        % split u and v
        U(:,n+1) = Y(1:length(u0),n+1);
        V(:,n+1) = Y(length(u0)+1:end, n+1);
    end
   
end
