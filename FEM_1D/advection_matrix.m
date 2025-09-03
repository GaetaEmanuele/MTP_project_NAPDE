function Matrix = advection_matrix(Dati,x0,x1,delta,T)
 
%=======================================================================================================
% This is the function where the matrix B is built, considering the fact
% that there are 2 types of tents (where the interval in space has lenght h
% and where the lenght is 2h)
%input: -Dati
%       - x0,x1 ending node
%       -delta = phi_top - pbi_bot
%       -Matrix classical FEM matrix 
%=======================================================================================================



    h = Dati.h;
    mu = Dati.mu;
    rho = Dati.rho;
    x_1 = x0:h:x1;
    l = length(x_1);
    
    if l == 3
        I1 = [x0,x_1(2)];
        I2 = [x_1(2),x0];
        tne = length(I1)-1;
        [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x_1(2));
    
        [basis] = shape_basis(Dati.element);

        [xq, W_1D] = Quadrature(nln);
        xq = sort(xq);
        nq = length(W_1D);
        [dphiq,Deriv] = evalshape(basis,xq);
        T1 = sparse(tnn,tnn);

        M1 = sparse(tnn,tnn); 
        for ie=1:tne
            iglo = connectivity(1:nln,ie);
            T_loc = zeros(nln,nln);

            M_loc = zeros(nln,nln);
            x_ = map(iglo,xq,x);
            for q=1:nq
                %B = BJ(q);
                for i=1:nln
                    x_q = x_(q);
                    phi_i = dphiq(i,q);
                    deriv_i = Deriv(1,i);
                    for j=1:nln
                        deriv_j = Deriv(1,j);
                        phi_j = dphiq(j,q);
                        T_loc(i,j) = T_loc(i,j) + delta(1,x_q)*deriv_i*phi_j*W_1D(q);
                        M_loc(i,j) = M_loc(i,j) + (h/2)*(phi_i*phi_j)*W_1D(q);
                        
                    end
                end
            end
                T1(iglo,iglo) = T1(iglo,iglo) + T_loc;

                M1(iglo,iglo) = M1(iglo,iglo) + M_loc;
        end
        
        tne = length(I2)-1;
        [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x_1(2), x1);
    
        [basis] = shape_basis(Dati.element);

        [xq, W_1D] = Quadrature(nln);
        xq = sort(xq);
        nq = length(W_1D);
        [dphiq,Deriv] = evalshape(basis,xq);
        T2 = sparse(tnn,tnn);

        M2 = sparse(tnn,tnn);
        for ie=1:tne
            iglo = connectivity(1:nln,ie);
            T_loc = zeros(nln,nln);

            M_loc = zeros(nln,nln);
            x_ = map(iglo,xq,x);
            for q=1:nq
                %B = BJ(q);
                for i=1:nln
                    phi_i = dphiq(i,q);
                    deriv_i = Deriv(1,i);
                    x_q = x_(q);
                    for j=1:nln
                        deriv_j = Deriv(1,j);
                        phi_j = dphiq(j,q);
                        T_loc(i,j) = T_loc(i,j) + delta(2,x_q)*deriv_i*phi_j*W_1D(q);
                        M_loc(i,j) = M_loc(i,j) + (h/2)*(phi_i*phi_j)*W_1D(q);
                        
                        
                    end
                end
            end
                T2(iglo,iglo) = T2(iglo,iglo) + T_loc;

                M2(iglo,iglo) = M2(iglo,iglo) + M_loc;
        end
        T(1:2, 1:2) = T(1:2, 1:2) + T1;
        T(2:3, 2:3) = T(2:3, 2:3) + T2;

        Matrix.T = T;

        Matrix.M1 = M1;
        Matrix.M2 = M2;
    else
        tne = length(x_1)-1;
        [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x1);
    
        [basis] = shape_basis(Dati.element);

        [xq, W_1D] = Quadrature(nln);
        xq = sort(xq);
        nq = length(W_1D);
        [dphiq,Deriv] = evalshape(basis,xq);
        T1 = sparse(tnn,tnn);
        for ie=1:tne
            iglo = connectivity(1:nln,ie);
            T_loc = zeros(nln,nln);

            x_ = map(iglo,xq,x);
            for q=1:nq
                %B = BJ(q);
                for i=1:nln
                    deriv_i = Deriv(1,i);
                    x_q = x_(q);
                    for j=1:nln
                        deriv_j = Deriv(1,j);
                        phi_j = dphiq(j,q);
                        T_loc(i,j) = T_loc(i,j) + delta(x_q)*deriv_i*phi_j*W_1D(q);
                        
                    end
                end
            end
                T(iglo,iglo) = T(iglo,iglo) + T_loc;
        end
    end
    Matrix.T = T;
end