function T = advection_matrix(Dati,x0,x1,delta,T)
    h = Dati.h;
    x_1 = x0:h:x1;
    l = length(x_1);
    if l == 3
        I1 = [x0,x(2)];
        I2 = [x(2),x0];
        tne = length(I1)-1;
        [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x(2));
    
        [basis] = shape_basis(Dati.element);

        [xq, W_1D] = Quadrature(nln);
        xq = sort(xq);
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
                    x_q = x(iglo(q));
                    for j=1:nln
                        phi_j = dphiq(j,q);
                        T_loc(i,j) = T_loc(i,j) + delta(x_q)*deriv_j*phi_j*W_1D(q);
                    end
                end
            end
                T1(iglo,iglo) = T1(iglo,iglo) + T_loc;
        end
        
        tne = length(I2)-1;
        [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x(2), x1);
    
        [basis] = shape_basis(Dati.element);

        [xq, W_1D] = Quadrature(nln);
        xq = sort(xq);
        [dphiq,Deriv] = evalshape(basis,xq);
        T2 = sparse(tnn,tnn);
        for ie=1:tne
            iglo = connectivity(1:nln,ie);
            T_loc = zeros(nln,nln);
            x_ = map(iglo,xq,x);
            for q=1:nq
                %B = BJ(q);
                for i=1:nln
                    deriv_i = Deriv(1,i);
                    x_q = x(iglo(q));
                    for j=1:nln
                        phi_j = dphiq(j,q);
                        T_loc(i,j) = T_loc(i,j) + delta(x_q)*deriv_j*phi_j*W_1D(q);
                    end
                end
            end
                T2(iglo,iglo) = T2(iglo,iglo) + T_loc;
        end
        T(1:2, 1:2) = T(1:2, 1:2) + T1;
        T(2:3, 2:3) = T(2:3, 2:3) + T2;

    else
        tne = length(x_1)-1;
        [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x1);
    
        [basis] = shape_basis(Dati.element);

        [xq, W_1D] = Quadrature(nln);
        xq = sort(xq);
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
                    x_q = x(iglo(q));
                    for j=1:nln
                        phi_j = dphiq(j,q);
                        T_loc(i,j) = T_loc(i,j) + delta(x_q)*deriv_j*phi_j*W_1D(q);
                    end
                end
            end
                T(iglo,iglo) = T(iglo,iglo) + T_loc;
        end
    end
end