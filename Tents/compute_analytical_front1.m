function [front,L] = compute_analytical_front1(x,Dati,t)
%=======================================================================================================
% this is the function where the map from the physical tent to reference
% one is built
%=======================================================================================================
% input = - x spatial node
%         - Dati 
%         - t the advacing front of the tent
% ouptu = -front : is a strcuture containing as function handle
% delta,phi,phix,phit


    %load data
    T = Dati.T;
    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    h = Dati.h;

    %evaluation of angular coefficient
    L = 0;
    if sum(t == 0)==2 || sum((t-T)<1e-3)==2
        L = max(t) - min(t);
        
    else
        L = (max(t) - min(t))/2; 
    end

    m = L/h;
    if sum(x==x0)==2
        if sum(t==0)==2
            phi_top = @(x) -m*x + L;
            delta = @(x) phi_top(x);
            phix = @(t) -m*t;
            phixt = -m;
            ft = @(x,t) t * phi_top(x);
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        elseif sum(abs(t-T)<1e-3)==2
            phi_bot = @(x) m*x + t(1);
            delta = @(x) T - phi_bot(x);
            ft = @(x,tau) tau * T + (1-tau)*phi_bot(x);
            phix = @(t) m * (1-t);
            phixt = -m;
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        else
            phi_top = @(x) -m*x + t(3);
            phi_bot = @(x) m*x + t(1);
            delta = @(x) phi_top(x) - phi_bot(x);
            ft = @(x,tau) tau * phi_top(x) + (1-tau)*phi_bot(x);
            phix = @(t) -m*t + (1-t)*m;
            phixt = -2*m;
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        end
    elseif sum(x==x1)==2
        if sum(t==0)==2
            phi_top = @(x_) m*(x_ - x(1)) + t(1);
            delta = @(x) phi_top(x);
            ft = @(x,tau) tau*phi_top(x);
            phix = @(tau) tau*m;
            phixt = m;
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        elseif sum(abs(t-T)<1e-3)==2
            phi_bot = @(x_) -m*(x_ - x(1)) + t(1);
            delta = @(x) T - phi_bot(x);
            ft = @(x,tau) tau*T + (1-tau)*phi_bot(x);
            phix = @(tau) -m*(1-tau);
            phixt = m;
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        else
            phi_top = @(x_) m*(x_-x(3)) + t(3);
            phi_bot = @(x_) -m*(x_ - x(2)) + t(2);
            delta = @(x) phi_top(x) - phi_bot(x);
            ft = @(x,t) t*phi_top(x) + (1-t)*phi_bot(x);
            phix = @(t) t*m + (1-t)* (-m);
            phixt = 2*m;
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        end
    else
        if sum(t==0)==2
            phi_top = @(a,x_) (a==1)*(m*(x_ - x(1))+ t(1)) + (a==2)*(-m*(x_ - x(3)) + t(3));
            delta = @(a,x) phi_top(a,x);
            ft = @(x_,t) t * phi_top(1,x_).*(x_ <= x(3)) + t * phi_top(2,x_).*(x_ > x(3));
            phix = @(a,t)  (a==1)*t*m + (a==2)*t*(-m);
            phixt = @(a) (a==1)*m + -m * (a==2);
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        elseif sum(abs(t-T)<1e-3)==2
            phi_bot = @(a,x_) (a==1)*(-m * (x_ - x(1)) + t(1)) + (a==2)*(m*(x_ - x(2))+ t(2));
            delta = @(a,x) T - phi_bot(a,x);
            ft = @(x_,tau) tau * T + (1-tau).*(x_ <= x(2))*phi_bot(1,x_) + (1-tau).*(x_ > x(2))*phi_bot(2,x_);
            phix = @(a,t) (1-t)*(a==1)*(-m) + (1-t)*(a==2)*m;
            phixt = @(a) m*(a==1) + -m*(a==2);
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        else
            phi_top = @(a,x_) (a==1)*(m*(x_ - x(1))+t(1)) + (a==2)*(-m*(x_ - x(4))+ t(4));
            phi_bot = @(a,x_) (a==1)*(-m*(x_ - x(1))+t(1)) + (a==2)*(m*(x_ - x(2))+t(2));
            delta = @(a,x) phi_top(a,x) - phi_bot(a,x);
            ft = @(x_,tau) tau * phi_top(1,x_).*(x_ <= x(2)) + tau * phi_top(2,x_).*(x_ > x(2)) + (1-tau)*phi_bot(1,x_).*(x_ <= x(2))+(1-tau)*phi_bot(2,x_).*(x_ > x(2));
            phix = @(a,t) ( t * m + (1-t)*(-m))*(a==1) + (-m*t + (1-t)*m)*(a==2);
            phixt = @(a) (a==1)*2*m -2*m*(a==2);
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        end
    end