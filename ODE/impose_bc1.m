function Y = impose_bc1(x,tau,Dati,front)
    
%=======================================================================================================
%Here are evaluated the BC taking as input x (spatial node), tau mapped
%time, Dati, front (function handle of the front)
%=======================================================================================================


    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    
    ft = front.ft;
    delta = front.delta;
    phix = front.phix;


    dim = 2*length(x);
    mu = Dati.mu;
    rho = Dati.rho;
    [~,t0] = direct_map1(x0,tau,ft);
    [~,t1] = direct_map1(x1,tau,ft);

    bcu =  inline(Dati.derivativex,'x','t');
    bcv = inline(Dati.derivativet,'x','t'); 
    Y = zeros(dim,1);

    if sum(x==x1)>0
        Y(dim/2) =  bcu(x1,t1) ;
        Y(end) = bcv(x1,t1);

    elseif sum(x==x0)>0
        Y(1) =  bcu(x0,t0) ;
        Y(dim/2+1) = bcv(x0,t0);
    end

end