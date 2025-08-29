function [front] = compute_analytical_front(x,Dati,t)
    %load data
    mu = Dati.mu;
    rho = Dati.rho;
    c = sqrt(mu/rho);
    T = Dati.T;

    x0 = Dati.domain(1);
    x1 = Dati.domain(2);
    x_min = min(x);
    x_max = max(x);
    if sum(x==x0) ==2
        b = t(3)+1/c * x(3);
        phi_top = @(x) -1/c * x + b;
        if sum(t==0) ==2
            delta = @(x) phi_top(x);
            phix = @(t) t*(-1/c);
            phixt = -1/c;
            ft = @(x,t) t*phi_top(x);
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        elseif sum(t==T) ==2
             b1 = -1/c * x(1) + t(1);
             phi_bot = @(x) 1/c *x + b1;
             delta = @(x) T-phi_bot(x);
             phix = @(t) (1-t)*1/c;
             phixt = -1/c;
             ft = @(x,t) (t)*T + (1-t)*phi_bot(x);
             front.ft = ft;
             front.delta = delta;
             front.phix = phix;
             front.phixt = phixt;
        else
             b1 = -1/c * x(1) + t(1);
             phi_bot = @(x) 1/c *x + b1;
             delta = @(x) phi_top(x) - phi_bot(x);
             phix = @(t) (1-t)*(1/c) + -t * 1/c;
             phixt =  -2/c;
             ft = @(x,t) (t)*phi_top(x) + (1-t)*phi_bot(x);
             front.ft = ft;
             front.delta = delta;
             front.phix = phix;
             front.phixt = phixt;
        end
    elseif sum(x==x1) == 2
            b = -1/c * x(3) + t(3);
            b1 = 1/c * x(2) + t(2);
            phi_top = @(x) 1/c * x +b;
            phi_bot = @(x) -1/c*x + b1;
            ft = @(x,t) (t)*phi_top(x) + (1-t)*phi_bot(x);
            delta = @(x) phi_top(x) - phi_bot(x);
            phix = @(t) (t)*(1/c) + (1-t)*(-1/c);
            phixt = 2/c;
            front.ft = ft;
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
    else
        if sum(t==0)==2
            b1 = -1/c * x(1) + t(1);
            b2 = 1/c * x(2) + t(2);
            phi_top = @(a,x)(1/c * x + b1)*(a==1) + (-1/c * x + b2)*(a==2);
            delta = @(a,x) phi_top(a,x);
            ft = @(x,t) t*(1/c * x + b1).*(x<=x(3)) + t*(-1/c * x + b2).*(x>x(3));
            front.ft = ft;
            phix = @(a,t) (t)*(1/c)*(a==1) + ((t)*(-1/c))*(a==2);
            phixt = @(a) 1/c * (a==1) + -1/c * (a==2);
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        elseif sum(t==T) ==2
            b1 = +1/c * x(1) + t(1);
            b2 = -1/c * x(3) + t(3);
            phi_bot = @(a,x)(-1/c * x + b1)*(a==1) + (+1/c * x + b2)*(a==2);
            delta = @(a,x) T - phi_bot(a,x);
            ft = @(x,t) (t)*T + (1-t)*(-1/c * x + b1).*(x<=x(2)) + (1-t)*(+1/c * x + b2).*(x>x(2));
            front.ft = ft;
            phix = @(a,t) (1-t)*(-1/c)*(a==1) + ((1-t)*(+1/c))*(a==2);
            phixt = @(a) 1/c * (a==1) + -1/c * (a==2);
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        else
            b1 = -1/c * x(1) + t(1);
            b2 = +1/c * x(3) + t(3);
            b3 = +1/c * x(1) + t(1);
            b4 = -1/c * x(3) + t(3);
            phi_top = @(a,x)(1/c * x + b1)*(a==1) + (-1/c * x + b2)*(a==2);
            phi_bot = @(a,x)(-1/c * x + b3)*(a==1) + (+1/c * x + b4)*(a==2);
            ft = @(x,t) (t)*(1/c * x + b1).*(x<=x(2)) + (t)*(-1/c * x + b2).*(x>x(2)) + (1-t)*(-1/c * x + b3).*(x<=x(2)) + (1-t)*(+1/c * x + b4).*(x>x(2));
            front.ft = ft;
            delta = @(a,x) phi_top(a,x) - phi_bot(a,x);
            %phix = @(a,t)(t*(1/c)+ (1-t)*(-1/c))*(a==1) +  (-t*(1/c)+ (1-t)*(1/c))*(a==2);
            phix = @(a,t)((2*t - 1)/c)*(a==1) + ((1 - 2*t)/c)*(a==2);

            phixt = @(a) 2/c * (a==1) + -2/c * (a==2);
            front.delta = delta;
            front.phix = phix;
            front.phixt = phixt;
        end
                
    
    end

end