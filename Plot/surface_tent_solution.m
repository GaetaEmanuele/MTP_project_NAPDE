function surface_tent_solution(x, front, Uhist, Vhist, Dati)
    dtau = Dati.dt;
    tau = 0:dtau:1;
    
    x = unique(x);

    Nx = numel(x);
    Nt = numel(tau);

    % costruzione griglia (X = nodi spaziali ripetuti, T = tempi fisici mappati)
    X = repmat(x(:)', Nt, 1);    % Nt × Nx
    T = zeros(Nt, Nx);
    for i = 1:Nt
        T(i,:) = front.ft(x, tau(i));   % tempo fisico per ogni x
    end

    % grafico
    figure
    tiledlayout(1,2,"Padding","compact","TileSpacing","compact")

    nexttile
    surf(T,X, Uhist'); shading interp
    xlabel('x'); ylabel('t'); zlabel('u'); title('uh')

    nexttile
    surf(T,X, Vhist'); shading interp
    xlabel('x'); ylabel('t'); zlabel('v'); title('vh')
end
