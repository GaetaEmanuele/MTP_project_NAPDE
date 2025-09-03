function [Matrices] = compute_matrix1(Dati,Matrix,front,x,t)
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
    
    T = Matrix.T;

    if N >2
        M1 = Matrix.M1;
        M2 = Matrix.M2;

        M12 = zeros(3,3);
        M21 = zeros(3,3);

        M12(1:2,1:2) = mu/rho * phi_x(1,t)*M1;
        M12(2:3,2:3) = M12(2:3,2:3) + mu/rho * phi_x(2,t)*M2;

        M21(1:2,1:2) = phi_x(1,t)*M1;
        M21(2:3,2:3) = M21(2:3,2:3) +  phi_x(2,t)*M2;
        
        M(1:N,N+1:end) = M12;
        M(N+1:end,1:N) = M21;

        S12 = zeros(3,3);
        S21 = zeros(3,3);

        S12(1:2,1:2) = mu/rho * phi_x_t(1)*M1;
        S12(2:3,2:3) = S12(2:3,2:3) + mu/rho * phi_x_t(2)*M2;

        S21(1:2,1:2) = phi_x_t(1)*M1;
        S21(2:3,2:3) = S21(2:3,2:3) +  phi_x_t(2)*M2;

        Sigma(1:N,N+1:end) = S12 + mu/rho * T;
        Sigma(N+1:end,1:N) = S21 + T;

    else
        M12 = mu/rho * phi_x(t)*M_;
        M21 = phi_x(t)*M_;
        M(1:N,N+1:end) = M12;
        M(N+1:end,1:N) = M21;
        S12 = mu/rho * (phi_x_t*M_ + T);
        S21 = phi_x_t*M_ + T;
        Sigma(1:N,N+1:end) = S12;
        Sigma(N+1:end,1:N) = S21;
    end
    
   
    Matrices.M = M;
    Matrices.Sigma = Sigma;


end