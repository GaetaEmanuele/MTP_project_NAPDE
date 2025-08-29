clear all
close all
clc
Dati=C_dati();

addpath FEM_1D

[Matrices,x]=C_matrix(Dati,0,1);
A = Matrices.A;
b = Matrices.f;
% Imposta Dirichlet omogenee sul nodo iniziale
A(1, :) = 0;
A(:, 1) = 0;
A(1, 1) = 1;

% Imposta Dirichlet omogenee sul nodo finale
A(end, :) = 0;
A(:, end) = 0;
A(end, end) = 1;

% Modifica anche il vettore dei termini noti
f(1) = 0;
f(end) = 0;
u = A\b;
uex = inline(Dati.exact_sol,'x');
figure()
plot(x,u)
hold on
plot(x,uex(x))
legend('uh','uex')