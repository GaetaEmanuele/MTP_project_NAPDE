function plot_u_v(nodes,Dati,U,V,t,Neven)
    x = zeros(length(U),1);
    if length(U)==Neven
        x = nodes(1,1:2:end);
        
    else
        x = nodes(1,2:2:end-1);
    end
    
    u = [];
    v = [];
    for i =1:length(U)
        u = [u,U{i}];
        v = [v,V{i}];
    end
    
    figure()
    dv = inline(Dati.derivativet, 'x', 't'); % se è diversa da du
    du = inline(Dati.derivativex, 'x', 't');  
    x_ = 0;
    if length(x)==Neven
        x_ = linspace(Dati.domain(1),Dati.domain(2),101);
    else
        x_ = linspace(Dati.domain(1)+ Dati.h,Dati.domain(2)-Dati.h,101);
    end
    u_interp = interp1(x, u, x_, 'spline');
    v_interp = interp1(x, v, x_, 'spline');
    duex = du(x_,t);
    dvex = dv(x_,t);
    subplot(2,1,1)
    plot(x_,u_interp)
    hold on 
    plot(x_,duex,'--')
    subplot(2,1,2)
    plot(x_,v_interp)
    hold on 
    plot(x_,dvex,'--')
end