function [nq, W_1D] = Quadrature(n)
    % Calcola i punti e i pesi di quadratura di Gauss-Legendre per n punti
    % n: numero di punti di quadratura
    % points: punti di quadratura
    % weights: pesi di quadratura

    % Utilizza funzioni di MATLAB per calcolare i punti e i pesi
    [nq, W_1D] = lgwt(n, -1, 1);
end

