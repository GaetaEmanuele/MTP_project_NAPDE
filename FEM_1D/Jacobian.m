function B = Jacobian(x, nq)
    % Calcola il Jacobiano in modo generale per elementi P1 o P2
    % x: coordinate fisiche dei nodi dell'elemento
    % nq: numero di nodi di quadratura
    % Output: B (Jacobiano, vettore di dimensione nq per P2 o scalare per P1)

    nln = length(x);  % Numero di nodi dell'elemento

    if nln == 2
        % Elemento P1 (lineare), Jacobiano costante
        B_ = (x(2) - x(1)) / 2;
        B = [B_,B_];
    elseif nln == 3
        % Elemento P2 (quadratico), Jacobiano variabile
        B = zeros(nq, 1);  % Pre-alloca il Jacobiano
        for q = 1:nq
            % Calcola il Jacobiano per ciascun nodo di quadratura
            % Ad esempio, puoi interpolare tra i nodi
            B(q) = 0.5 * ((x(2) - x(1)) + (x(3) - x(1)) * (2 * q / nq - 1)); 
        end
    else
        error('Tipo di elemento non supportato');
    end
end

