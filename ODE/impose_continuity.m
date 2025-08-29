function [u0,v0] = impose_continuity(dependence,list_U,list_V,x_hat,Dati)
    %select the tents
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    index_tent = dependence;
    index_tent = index_tent(index_tent>0);
    U_selected = list_U(index_tent);
    V_selected = list_V(index_tent);
    u0 = zeros(1,length(x_hat));
    v0 = zeros(1,length(x_hat));
    %from x extract the medium indexes
    
    index = 2;
    %distiction between triangular and quadrilateral tents
    if length(U_selected)==2
        %aggragate the initial condition from the 2 tents
        U1 = U_selected{1};
        U2 = U_selected{2};
        V1 = V_selected{1};
        V2 = V_selected{2};
        if (size(U1,1)==2) 
            u0 = [U1(:,end);U2(2:index)];
            v0 = [U1(:,end);U2(2:index)];
        else 
            u0 = [U1(index:end,end);U2(2:index,end)];
            v0 = [V1(index:end,end);V2(2:index,end)];
        end
    else
        U = U_selected{1};
        V = V_selected{1};
        if min(x_hat) == x0
            u0 = U(1:index,end);
            v0 = V(1:index,end);
        elseif max(x_hat) == x1
            u0 = U(index:end,end);
            v0 = V(index:end,end);
        end
    end
end