function [Ui,Vi] = solver(Dati,dependence,U,V,x_hat,t_hat,front)
%discritization in space
[Matrices,x]=C_matrix(Dati,min(x_hat),max(x_hat),front.delta);

%solve ODE
[Ui,Vi] = runge_kutta_wave_eq(Dati,x,dependence,Matrices,U ,V,front,t_hat);
%[Ui,Vi] = BE_wave_eq(Dati,x,dependence,Matrices,U ,V,front,t_hat);
%[Ui,Vi] = LDDRK_wave_eq(Dati,x,dependence,Matrices,U ,V,front,t_hat);
%[Ui,Vi] = runge_kutta38_wave_eq(Dati,x,dependence,Matrices,U ,V,front,t_hat);
end
