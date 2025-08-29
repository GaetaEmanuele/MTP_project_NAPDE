function [u,v] = map_solution(U_hat_i,V_hat_i,front,x,t,Dati)
delta = front.delta;
phix = front.phix;
n = size(U_hat_i,1);
x0 = Dati.domain(1);
xN = Dati.domain(2);

    if n==2
        if sum(abs(x-x0)<1e-3)>0
            u = U_hat_i(1,end);
            v = V_hat_i(1,end);

        else
            u= U_hat_i(2,end);
            v = V_hat_i(2,end);

        end
    else
        if length(x)==3
            x_ = x(3);
            u= U_hat_i(2,end);
            v = V_hat_i(2,end);
        else
            x_ = x(2);
            u = U_hat_i(2,end);
            v = V_hat_i(2,end);
        end

    end



end