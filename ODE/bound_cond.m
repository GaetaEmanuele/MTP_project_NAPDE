function [u_bd,v_bd] = bound_cond(Dati,x_bd,t,u_bd,v_bd)
    mu = Dati.mu;
    rho = Dati.rho;
    Funu = inline(Dati.derivativex,'x','t');
    Funv = inline(Dati.derivativet,'x','t');
    if sum(x_bd==Dati.domain(1))>0
        u_bd(1) = - mu * Funu(Dati.domain(1),t);
        v_bd(1)= rho * Funv(Dati.domain(1),t);
    elseif sum(x_bd==Dati.domain(2))>0
        u_bd(end) = - mu *Funu(Dati.domain(2),t);
        v_bd(end)= rho * Funv(Dati.domain(2),t);
    end

end