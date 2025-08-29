function [Matrices] = compute_matrix(Dati,Matrix,front,x,t)
    %load the data
    mu = Dati.mu;
    rho = Dati.rho;
    x0 = Dati.domain(1);
    xN= Dati.domain(2);
    %load as function handle dphi/dx and d(dphi/dx)/dt
    phi_x = front.phix;
    phi_x_t = front.phixt;


    M_ = Matrix.M;
    N = size(M_,1);
    M = zeros(2*N,2*N);
    Sigma = zeros(2*N,2*N);
    M(1:N,1:N) = M_;
    M(N+1:end,N+1:end) = M_;

    if N>2
        T1 = + mu/rho * Matrix.M1 * phi_x(1,t);
        T2 = + mu/rho * phi_x(2,t)* Matrix.M2;
        A12 = zeros(N,N);
        A12(1:2,1:2)= A12(1:2,1:2)+T1;
        A12(2:3,2:3) = A12(2:3,2:3)+T2;
        M(1:N,N+1:end) =  A12;
        T1 = + Matrix.M1 * phi_x(1,t);
        T2 = + phi_x(2,t)* Matrix.M2;
        A21 = zeros(N,N);
        A21(1:2,1:2)= A21(1:2,1:2)+T1;
        A21(2:3,2:3) = A21(2:3,2:3)+T2;
        M(N+1:end,1:N) =  A21;
        
        

        Tvol1 = mu/rho * Matrix.T;
        Tvol2 = Matrix.T;
        T1 = mu/rho * Matrix.M1 * phi_x_t(1);
        T2 = mu/rho * phi_x_t(2)* Matrix.M2;
        Tvol1(1:2,1:2) = Tvol1(1:2,1:2)+T1;
        Tvol1(2:3,2:3) = Tvol1(2:3,2:3)+T2;
        
        T1 = Matrix.M2 * phi_x_t(1);
        T2 = phi_x_t(2)* Matrix.M2;
        Tvol2(1:2,1:2) = Tvol2(1:2,1:2)+T1;
        Tvol2(2:3,2:3) = Tvol2(2:3,2:3)+T2;
        Sigma(1:N,N+1:end) = Tvol1;
        Sigma(N+1:end,1:N) = Tvol2;

    else
        T1 =  +mu/rho * phi_x(t)*Matrix.M;
        A12 = T1;
        M(1:N,N+1:end) = A12;
        T1 = Matrix.M * phi_x(t);
        A21 = zeros(N,N);
        A21= T1;
        M(N+1:end,1:N) = A21;
        Tvol1 = mu/rho * Matrix.T;
        Tvol2 = Matrix.T;

        T1 = mu/rho * Matrix.M * phi_x_t;
        Tvol1 = Tvol1+T1;
        
        T1 = Matrix.M * phi_x_t;
        Tvol2=Tvol2 + T1;
        Sigma(1:N,N+1:end) = Tvol1;
        Sigma(N+1:end,1:N) = Tvol2;
    end

     %impose that the solution is constant on the chart. lines
    l = zeros(2*N,1);
    lt = l';
    if sum(x == x0)>0
        M(:,1) = l;
        M(1,:) = lt;
        M(:,N+1) = l;
        M(N+1,:) = lt;
        M(1,1) = 1;
        M(N+1,N+1) = 1;

        Sigma(:,1) = l;
        Sigma(1,:) = lt;
        Sigma(:,N+1) = l;
        Sigma(N+1,:) = lt;
        Sigma(1,N+1) = 1;
        Sigma(N+1,1) = 1;
    elseif sum(x==xN)>0
        M(:,N) = l;
        M(N,:) = lt;
        M(:,end) = l;
        M(end,:) = lt;
        M(N,N) = 1;
        M(end,end) = 1;

        Sigma(:,N) = l;
        Sigma(N,:) = lt;
        Sigma(:,end) = l;
        Sigma(end,:) = lt;
        Sigma(N,end) = 1;
        Sigma(end,N) = 1;
    end

   

    Matrices.M = M;
    Matrices.Sigma = Sigma;

end