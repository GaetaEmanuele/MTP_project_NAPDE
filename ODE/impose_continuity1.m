function [u0,v0] = impose_continuity1(dependence,list_U,list_V,x_hat,Dati)
    tol = 1e-4;
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    x_max = max(x_hat);
    x_min = min(x_hat);
    h = Dati.h;
    index_tent = dependence;
    index_tent = index_tent(index_tent>0);
    U_selected = list_U(index_tent);
    V_selected = list_V(index_tent);
    u0 = zeros(1,length(x_hat));
    v0 = zeros(1,length(x_hat));
    if length(x_hat)==2
        U = U_selected{1};
        V = V_selected{1};
        if sum(x_hat==x0)>0
            u0 = U(1:2,end);
            v0 = V(1:2,end);
        else
            u0 = U(2:3,end);
            v0 = V(2:3,end);
        end
    else
        U1 = U_selected{1};
        U2 = U_selected{2};
        V1 = V_selected{1};
        V2 = V_selected{2};
        if size(U1,1)==2
            u0 = [U1(:,end);U2(2,end)];
            v0 = [V1(:,end);V2(2,end)];
        else
            u0 = [U1(2:3,end);U2(2,end)];
            v0 = [V1(2:3,end);V2(2,end)]; 
        end

    end

end