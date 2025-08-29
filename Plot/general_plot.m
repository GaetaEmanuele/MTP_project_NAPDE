function general_plot(n_even,n_odd,U,V,Dati,nodes,delta)
    it = 1;
    N = n_even; 

    t = 0;
    T = Dati.T;

    while it < length(U)
        if t +delta < T
            t = t + delta;
        else
            t = T;
        end

        u = [];
        v = [];
        for i = it : it + N - 1
            u = [u, U{i}];
            v = [v, V{i}];
        end

        if N == n_even
            x = nodes(1,1:2:end);
            N = n_odd;
        else
            x = nodes(1,2:2:end-1);
            N = n_even;
        end

        it = i + 1;

        dv = inline(Dati.derivativet, 'x', 't'); % se diversa da du
        du = inline(Dati.derivativex, 'x', 't');

        if length(x) == n_even
            x_ = linspace(Dati.domain(1), Dati.domain(2), 101);
        else
            x_ = linspace(Dati.domain(1) + Dati.h, Dati.domain(2) - Dati.h, 101);
        end

        u_interp = interp1(x, u, x_, 'spline');
        v_interp = interp1(x, v, x_, 'spline');
        duex = du(x_, t);
        dvex = dv(x_, t);

        clf; % Pulisce la figura ogni volta
        subplot(2,1,1)
        plot(x_, u_interp, 'b', 'DisplayName', 'u numerica')
        hold on
        plot(x_, duex, 'r--', 'DisplayName', 'u esatta')
        title(['Tempo t = ', num2str(t)])
        legend show
        grid on

        subplot(2,1,2)
        plot(x_, v_interp, 'b', 'DisplayName', 'v numerica')
        hold on
        plot(x_, dvex, 'r--', 'DisplayName', 'v esatta')
        legend show
        grid on

        drawnow;
        pause(0.05) % rallenta l'animazione
    end
end
