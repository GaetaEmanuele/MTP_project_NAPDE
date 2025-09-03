function [Matrices,x]=C_matrix(Dati,x_0,x_1,delta)
    
%=======================================================================================================
% This is the function where the Mass matrix and the matrix B are built
%input: -Dati
%       - x_0,x_1 the ending node
%       -delta = phi_top - phi_bot 
%outpu: -Matrices (strcuture containing the matrix M,B)
%=======================================================================================================



    x0 = x_0;
    x1 = x_1;
    h = Dati.h;
    x_temporary = x0:h:x1;
    tne = length(x_temporary)-1;
    %creating mesh
    [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x1);
    
    %basis function
    [basis] = shape_basis(Dati.element);
    
    %quadrature point
    [xq, W_1D] = Quadrature(nln);
    xq = sort(xq);
    %evaluation basis function
    [dphiq,Deriv] = evalshape(basis,xq);

    
    M = sparse(tnn,tnn);
    
    nq = length(W_1D);
    %BJ = Jacobian(xq,nq);

    % Initialization of Matrices
    Matrices.M = M;

    for ie=1:tne
        iglo = connectivity(1:nln,ie);

        
        M_loc = zeros(nln,nln);
        
        
        x_ = map(iglo,xq,x);
        for q=1:nq
            %B = BJ(q);
            for i=1:nln
                deriv_i = Deriv(1,i);
                phi_i = dphiq(i,q);
                x_q = x(iglo(q));
                
                for j=1:nln
                    deriv_j = Deriv(1,j);
                    phi_j = dphiq(j,q);
                    M_loc(i,j) = M_loc(i,j) + (h/2)*(phi_i*phi_j)*W_1D(q);
                end
            end
        end
        M(iglo,iglo) = M(iglo,iglo) + M_loc;
    end
    % store the matrices
    Matrices.M = M;
    T = sparse(tnn,tnn);
    Matrix = advection_matrix(Dati,x0,x1,delta,T);
    Matrices.T = Matrix.T;
    if length(x) == 3
        Matrices.M1 = Matrix.M1;
        Matrices.M2 = Matrix.M2;
    end
end