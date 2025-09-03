function inflow = numerical_flux(front,Dati,x)
    
%=======================================================================================================
% Here is built the numerical Flux term
%input: - front strcuture of function handle that describe the map
%       - Dati
%       - x spatial coordiante
%=======================================================================================================
%

    delta = front.delta;
    phixt = front.phixt;
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    dim = 2*length(x);
    mu = Dati.mu;
    rho = Dati.rho;
    inflow = zeros(dim,dim);
    if dim==6
        inflow(1,dim/2 +1) = -  delta(1,min(x))* mu/rho;
        inflow(dim/2 +1,1) = - delta(1,min(x));
        
        inflow(dim/2,end) = delta(2,max(x))* mu/rho;
        inflow(end,dim/2) = delta(2,max(x));
    else
        inflow(1,dim/2 +1) = -  delta(min(x))* mu/rho;
        inflow(dim/2 +1,1) = - delta(min(x));
        
        inflow(dim/2,end) = delta(max(x))* mu/rho;
        inflow(end,dim/2) = delta(max(x));
       
    end

end