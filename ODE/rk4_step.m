function k = rk4_step(Y, t,x,Dati,Matrix,front)
        % evaluate the time dependent matrixes
        Matrices = compute_matrix1(Dati,Matrix,front,x,t);
        M = Matrices.M;
        W = Matrices.W;
        Sigma = Matrices.Sigma;
        D = Matrices.D;
        
        eps = 0.5*1e-2*Dati.h*sqrt((D*Y)'*W*(D*Y));
        eps_min = 1e-6;
        eps_max = 1e-2;
        eps = max(min(eps, eps_max), eps_min);

        
        dim = length(Y);
        I = numerical_flux(front,Dati,x);

        % Derivate
        dy_dt = (M)\(-Sigma*Y + I*Y);
        
        
        
        % Concatenation
        k = dy_dt;
end