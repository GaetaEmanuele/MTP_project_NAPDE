function Y = impose_bc(x,tau,Dati,front)
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
        if dim==4
            Y(dim/2) =  mu *bcu(x1,t1) + mu * bcv(x1,t1)*phix(tau);
            Y(end) = rho*delta(x1)*bcv(x1,t1);
        else
            Y(dim/2) =  mu *bcu(x1,t1) + mu * bcv(x1,t1)*phix(2,tau);
            Y(end) = rho*delta(2,x1)*bcv(x1,t1);
        end


    elseif sum(x==x0)>0
        if dim==4
            Y(1) =  mu *bcu(x0,t0) + mu * bcv(x0,t0)*phix(tau);
            Y(dim/2+1) = delta(x0)*bcv(x0,t0);
        else
            Y(1) =  mu *bcu(x0,t0) + mu * bcv(x0,t0)*phix(1,tau);
            Y(dim/2+1) = delta(1,x0)*bcv(x0,t0);
        end
    end

end