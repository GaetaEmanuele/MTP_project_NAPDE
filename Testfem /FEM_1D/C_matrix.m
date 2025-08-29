function [Matrices,x]=C_matrix(Dati,x_0,x_1,delta)
    x0 = x_0;
    x1 = x_1;
    h = Dati.h;
    x_temporary = x0:h:x1;
    tne = length(x_temporary)-1;
    [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x1);
    
    [basis] = shape_basis(Dati.element);

    [xq, W_1D] = Quadrature(nln);
    xq = sort(xq);
    [dphiq,Deriv] = evalshape(basis,xq);

    %A = sparse(tnn,tnn);
    M = sparse(tnn,tnn);
    f = sparse(tnn,1);
    nq = length(W_1D);
    %BJ = Jacobian(xq,nq);

    % Initialization of Matrices
    %Matrices.A = A;
    Matrices.M = M;
    Matrices.f = f;

    for ie=1:tne
        iglo = connectivity(1:nln,ie);

        %K_loc = zeros(nln,nln);
        M_loc = zeros(nln,nln);
        f_loc = zeros(nln,1);
        fun = inline(Dati.force,'x');
        x_ = map(iglo,xq,x);
        for q=1:nq
            %B = BJ(q);
            for i=1:nln
                deriv_i = Deriv(1,i);
                phi_i = dphiq(i,q);
                x_q = x(iglo(q));
                f_loc(i) = f_loc(i) + (h/2)*fun(x_q)*phi_i*W_1D(q);
                for j=1:nln
                    deriv_j = Deriv(1,j);
                    phi_j = dphiq(j,q);
                    %K_loc(i,j) = K_loc(i,j) + (2/h)*(deriv_j*deriv_i)*W_1D(q);
                    M_loc(i,j) = M_loc(i,j) + (h/2)*(phi_i*phi_j)*W_1D(q);
                end
            end
        end
        %A(iglo,iglo) = A(iglo,iglo) + K_loc ;
        M(iglo,iglo) = M(iglo,iglo) + M_loc;
        f(iglo) = f(iglo) + f_loc;
    end
    % store the matrices
    %Matrices.A = A;
    Matrices.M = M;
    Matrices.f = f;
    T = sparse(nln,nln);
    T = advection_matrix(Dati,x0,x1,delta,T);
    Matrices.T = T;
end