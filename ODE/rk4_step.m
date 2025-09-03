function k = rk4_step(Y, t,x,Dati,Matrix,front)
        
%=======================================================================================================
% This is the function where the Ki coefficent of the RK method are
% evaluated
%=======================================================================================================
%

        % evaluate the time dependent matrixes
        Matrices = compute_matrix1(Dati,Matrix,front,x,t);
        M = Matrices.M;

        Sigma = Matrices.Sigma;
        

        
        dim = length(Y);
        I = numerical_flux(front,Dati,x);

        % Derivate
        dy_dt = (M)\(-Sigma*Y + I*Y);
        
        
        
        % Concatenation
        k = dy_dt;
end