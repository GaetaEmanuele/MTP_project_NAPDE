function [u0,v0] = impose_ic(Dati,front,x)

%=======================================================================================================
% Here are evaluated the IC for the 1st layer of tent 
%input = Dati,front,x
%output = u0 itial uh, v0 itial vh
%=======================================================================================================
%

ft = front.ft;
phix = front.phix;
delta = front.delta;
mu = Dati.mu;
rho = Dati.rho;
x0 = Dati.domain(1);
x1 = Dati.domain(2);



u = inline(Dati.derivativex,'x','t');

v = inline(Dati.derivativet,'x','t');
[~,tau] = direct_map1(x,0,ft);
u0 =  zeros(length(x),1);
v0 =  zeros(length(x),1);

u0 = mu*u(x,tau);
 

if length(x)==3
    v0(1:2) = rho * delta(1,x(1:2)).*v(x(1:2),tau(1:2));
    v0(3) = rho * delta(2,x(3)).*v(x(3),tau(3));
else
    v0 = rho * delta(x).*v(x,tau);
end

end