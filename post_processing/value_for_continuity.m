function [u,v] = value_for_continuity(U,V,x,t,up_front,lower_front,delta,Dati)
    h = Dati.h;
    dt = Dati.dt;
    [x_hat,t_hat] = inverse_map(x,t,up_front,lower_front,delta);
    index_x = floor((x_hat-min(x_hat))/h + 1);
    index_t = floor(t_hat/dt + 1);
    u = zeros(length(index_x),1);
    v = zeros(length(index_x),1);
    for i=1 : length(index_x)
        u(i) = U(index_x(i),index_t(i));
        v(i) = V(index_x(i),index_t(i));
    end
end