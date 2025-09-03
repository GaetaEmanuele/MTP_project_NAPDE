function e = compute_error1(Dati, U, V, time)
%---------------------------------------------------------
% evaluation of the L2 combined error of u and v
%
% INPUT:
%   Dati : struct 
%   U    : nodal value of numerical solution u_h
%   V    : nodal value of numerical solution v_h
%   time : final time T
%
% OUTPUT:
%   e    : global error sqrt( ||u-uh||^2 + ||v-vh||^2 )
%---------------------------------------------------------
    
    t = time;
    h = Dati.h;
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    x_temporary = x0:h:x1;
    u_ex = inline(Dati.derivativex,'x','t');
    v_ex = inline(Dati.derivativet,'x','t');

    tne = length(x_temporary)-1;

    [ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x1);
   
   

    % quadratura di Gauss su [-1,1]
    [xi, w] = Quadrature(nln);

    e2_u = 0; 
    e2_v = 0;

    % loop over the elment
    for K = 1:tne
        xL = x(K); 
        xR = x(K+1); 
        hK = xR - xL;

        % nodal value of the numerical sol
        UhL = U(K);   UhR = U(K+1);
        VhL = V(K);   VhR = V(K+1);

        % loop over quadratur point
        for q = 1:length(w)
            % map over [xL, xR]
            xq = ((xR - xL)/2)*xi(q) + (xR + xL)/2;

            % local shape function
            phiL = (xR - xq)/hK;
            phiR = (xq - xL)/hK;

            % approx FEM u_h, v_h
            uhq = UhL*phiL + UhR*phiR;
            vhq = VhL*phiL + VhR*phiR;

            % exact sol
            uq = u_ex(xq, t);
            vq = v_ex(xq, t);

            % building the error
            e2_u = e2_u + (hK/2)*w(q) * (uq - uhq)^2;
            e2_v = e2_v + (hK/2)*w(q) * (vq - vhq)^2;
        end
    end

    % total error
    e = sqrt(e2_u + e2_v);

end
