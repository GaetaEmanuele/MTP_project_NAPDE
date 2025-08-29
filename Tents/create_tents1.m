function [tents,nodes,tentDependencies,n_even,n_odd] = create_tents1(Dati)
    %load the data

    x0 = Dati.domain(1);
    xN = Dati.domain(2);
    h = Dati.h;
    T = Dati.T;
    x = x0:h:xN;
    N = length(x);
    c = sqrt(Dati.mu / Dati.rho);
    %generate the advancing front 
    [nodes,n_even,n_odd] = C_nodes(x,T,c,N);
    %initialize tents 
    tents = [];
    %Aggregate the tents
    iter = 0;
    tol = 1e-12;  % se vuoi gestire approssimazioni numeriche

    % Cerco dalla seconda riga in poi (quindi righe 2:end)
    idx = find(all(abs(nodes(2:end,:) - T) < tol, 2), 1);
    Tent = zeros(4,2);
    Ncheck = n_even;
    for i=2:idx
        if Ncheck==n_even
            if iter==0
                iter = 1;
                x_cord = [x(1),x(2),x(1)];
                t = [nodes(i,1),nodes(i,2),nodes(i+1,1)];
                Tent = StructTent(x_cord, t);
                tents(:,:,1) = Tent;
                for j=2:2:(N-3)
                    x_cord = [x(j),x(j+2),x(j+1)];
                    t = [nodes(i,j),nodes(i,j+2),nodes(i+1,j+1)];
                    Tent = StructTent(x_cord, t);
                    tents(:,:,end+1) = Tent;
                end
                x_cord = [x(end-1),x(end),x(end)];
                t = [nodes(i,end-1),nodes(i,end),nodes(i+1,end)];
                Tent = StructTent(x_cord, t);
                tents(:,:,end+1) = Tent;
            else
                x_cord = [x(1),x(2),x(1)];
                t = [nodes(i-1,1),nodes(i,2),nodes(i+1,1)];
                Tent = StructTent(x_cord, t);
                tents(:,:,end+1) = Tent;
                for j=2:2:(N-3)
                    x_cord = [x(j),x(j+1),x(j+2),x(j+1)];
                    t = [nodes(i,j),nodes(i-1,j+1),nodes(i,j+2),nodes(i+1,j+1)];
                    Tent = StructTent(x_cord, t);
                    tents(:,:,end+1) = Tent;
                end
                x_cord = [x(end-1),x(end),x(end)];
                t = [nodes(i,end-1),nodes(i-1,end),nodes(i+1,end)];
                Tent = StructTent(x_cord, t);
                tents(:,:,end+1) = Tent;

            end
            Ncheck = n_odd;
        else
            for j=1:2:(N-2)
                    x_cord = [x(j),x(j+1),x(j+2),x(j+1)];
                    t = [nodes(i,j),nodes(i-1,j+1),nodes(i,j+2),nodes(i+1,j+1)];
                    Tent = StructTent(x_cord, t);
                    tents(:,:,end+1) = Tent;
            end
            Ncheck = n_even;
        end
     
    end

    %Creation of the tree of dependencies
    tentDependencies = findTentDependencies1(tents,n_even,n_odd);
end