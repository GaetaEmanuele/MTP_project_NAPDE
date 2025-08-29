function [U] = direct_map(x_hat,t_hat,U_hat_i,x0,x1,x,t)
   
    if length(U_hat_i)==2
        if min(x_hat)==x0
            U_i = interp1(t_hat(2:3), U_hat_i, t(2:3));
        elseif max(x_hat)==x1
            t_hat = [t_hat(1);t_hat(3)];
            t = [t(1);t(3)];
            U_i = interp1(t_hat, U_hat_i, t);
        end
        U = U_i(end);
    else 
        if length(x_hat)==3 
            F = scatteredInterpolant(x_hat, t_hat, U_hat_i);
            U_i = F(x,t);
        else 
            x_hat = [x_hat(1);x_hat(3);x_hat(4)];
            t_hat = [t_hat(1);t_hat(3);t_hat(4)];
            x = [x(1);x(3);x(4)];
            t = [t(1);t(3);t(4)];
            F = scatteredInterpolant(x_hat, t_hat, U_hat_i);
            U_i = F(x,t);
        end
        U = U_i(2);
    end
end