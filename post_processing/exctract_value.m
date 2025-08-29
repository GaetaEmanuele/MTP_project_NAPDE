function [value,x,t] = exctract_value(U_,x,t,up_front,lower_front,delta,Dati)
    tol = 1e-4;
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    h = Dati.h;
    dt = Dati.dt;
    if sum(t==0)>1
        if abs(x(1)-x0)>tol
            x = [x(1);x(end);x(2:end)];
            t = [t(1);t(1);t(2:end)];
        end
    elseif sum(abs(t-Dati.T)<tol)>1
        if abs(x(1)-x0)>tol
            t = [t(1);t(end);t(1)];
        else
            x = x(2:end);
            t = t(2:end);
        end
    else
        if abs(min(x)-x0)<tol && abs(max(x)-min(x)-h)<tol
            t = [t(2);t(2:end)];
        elseif abs(max(x)-x1)<tol && abs(max(x)-min(x)-h)<tol
            t = [t(1);t(1);t(end)];
        else
            t = [t(1);t(1);t(3:end)];
        end
    end
    [x_hat,t_hat] = inverse_map(x,t,up_front,lower_front,delta);
    index_x = floor((x_hat-min(x_hat))/h + 1);
    index_t = floor(t_hat/dt + 1);
    value = zeros(length(index_x),1);
    for i=1 : length(index_x)
        value(i) = U_(index_x(i),index_t(i));
    end
end