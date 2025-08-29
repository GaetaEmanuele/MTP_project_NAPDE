function e = compute_error1(Dati, U, V, time)
%---------------------------------------------------------
% Calcola l'errore L2 combinato di u e v al tempo finale
%
% INPUT:
%   Dati : struct con campi
%          - nodes      : vettore nodi globali (1 x (Ne+1))
%          - Ne         : numero di elementi
%          - nq         : numero punti quadratura per elemento
%          - u_exact(x,t)
%          - v_exact(x,t)
%   U    : valori nodali della soluzione numerica u_h
%   V    : valori nodali della soluzione numerica v_h
%   time : tempo finale T
%
% OUTPUT:
%   e    : errore totale sqrt( ||u-uh||^2 + ||v-vh||^2 )
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

    % loop sugli elementi
    for K = 1:tne
        xL = x(K); 
        xR = x(K+1); 
        hK = xR - xL;

        % valori nodali locali delle soluzioni numeriche
        UhL = U(K);   UhR = U(K+1);
        VhL = V(K);   VhR = V(K+1);

        % loop punti di quadratura
        for q = 1:length(w)
            % mappa punto di quadratura su [xL, xR]
            xq = ((xR - xL)/2)*xi(q) + (xR + xL)/2;

            % funzioni di forma P1 locali
            phiL = (xR - xq)/hK;
            phiR = (xq - xL)/hK;

            % approx FEM u_h, v_h
            uhq = UhL*phiL + UhR*phiR;
            vhq = VhL*phiL + VhR*phiR;

            % soluzioni esatte
            uq = u_ex(xq, t);
            vq = v_ex(xq, t);

            % contributi all'errore
            e2_u = e2_u + (hK/2)*w(q) * (uq - uhq)^2;
            e2_v = e2_v + (hK/2)*w(q) * (vq - vhq)^2;
        end
    end

    % errore totale
    e = sqrt(e2_u + e2_v);

end
