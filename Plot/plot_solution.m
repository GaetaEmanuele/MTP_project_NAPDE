function plot_solution(acceptable_x,t,U_hat_i,V_hat_i,Dati,front)
    
    % function where you can plot the solution inside the ith tent

    %plot of the velocity 
    h = Dati.h;
    x_ = min(acceptable_x):h:max(acceptable_x);

    % Definizione delle derivate (usa funzioni diverse se necessario)
    dv = inline(Dati.derivativet, 'x', 't'); % se è diversa da du
    du = inline(Dati.derivativex, 'x', 't');

for j = 1:length(t)
    tj = t(j);

    % Mappa inversa (o diretta)
    [x_phys, t_phys] = direct_map1(x_, tj, front.ft);

    % Valutazione delle derivate esatte
    dvex = dv(x_phys, t_phys);
    duex = du(x_phys, t_phys);

    % Plot per V
    figure(1)
    subplot(2,1,1)
    plot(x_, V_hat_i(:,j), 'r-', 'LineWidth', 1.5)
    ylabel('v numerica')
    title(['Tempo t = ', num2str(tj)])
    ylim([-pi, pi])

    subplot(2,1,2)
    plot(x_, dvex, 'b--', x_, V_hat_i(:,j), 'r-', 'LineWidth', 1.5)
    ylabel('v esatta')
    xlabel('x')
    ylim([-pi, pi])

    % Plot per U
    figure(2)
    subplot(2,1,1)
    plot(x_, U_hat_i(:,j), 'r-', 'LineWidth', 1.5)
    ylabel('u numerica')
    title(['Tempo t = ', num2str(tj)])
    ylim([-pi, pi])

    subplot(2,1,2)
    plot(x_, duex, 'b--', x_, U_hat_i(:,j), 'r-', 'LineWidth', 1.5)
    ylabel('u esatta')
    xlabel('x')
    ylim([-pi, pi])

    pause(0.05)
end

end

