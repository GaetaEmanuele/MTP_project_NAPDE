function [u] = Sommerfeld_bd_cond(Y,c)
    u = Y(1:length(u0), n+1);
    u(1) = u(2) - c * (u(2) - u(1)) / dx * dt; % Bordo sinistro
    u(end) = u(end-1) - c * (u(end) - u(end-1)) / dx * dt; % Bordo destro
end