function [matrices] = impose_sm_bd(matrices,bd_cond,Dati)
    A = matrices.A;
    M = matrices.M;
    n = size(A,1);
    D = zeros(n,n);
    if strcmp(bd_cond.field1,'Sommerfeld')
        A(1,1) = A(1,1) + Dati.c;
        M(1,1) = M(1,1) + Dati.c;
        D(1,1) = Dati.dumpling_factor;
    elseif strcmp(bd_cond.field2,'Sommerfeld')
        A(end,end) = A(end,end) + Dati.c;
        M(end,end) = M(end,end) + Dati.c;
        D(end,end) = Dati.dumpling_factor;
    end
    if strcmp(bd_cond.field1,'Sommerfeld') || strcmp(bd_cond.field2,'Sommerfeld')
        matrices.A = A;
        matrices.M = M;
        matrices.D = D;
    end

end