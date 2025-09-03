function [sol_u,sol_v,nodes,n_even,error] = c_main(Dati)

%=======================================================================================================
% This is the main file of the project where all the overall structure is
% conrdinated
%=======================================================================================================
% input = Dati the data structure containing the information necessary to
% solve the problem
% output = 
%         - sol_u is the numerical flux evaluated on the free node
%         - sol_v is the numerical velocity evaluated on the free node
%         - nodes is the structure containing the space-time node necessary for create tents         
%         - n_even 
%         - error is the following error e^2 = normL2 (u-uh) ^2 + normL2
%         (v-vh)^2


addpath Tents
addpath Plot
addpath map
addpath Dependeces
addpath ODE
addpath FEM_1D
addpath Solver
addpath post_processing
addpath Error

if nargin <1
    Dati  = C_dati();
end
fprintf('Creating tents...\n');
%creation of the space-time meshes
[tents,nodes,tentDependencies,n_even,n_odd] = create_tents1(Dati);



if Dati.visual_graph == 'Y'
    C_plot(nodes,Dati);
end
I_max = size(tents,3);
n_sol = size(nodes,1)-3;
n_col = size(nodes,2);
U = {};
V = {};
sol_u = {};
sol_v = {};
x0 = Dati.domain(1);
x1 = Dati.domain(2);
res = zeros(I_max,1);
N = n_even;
it =1;
time1 = tents(:,2,1);
delta = max(time1);
list_u = [];
list_v = [];
T =0;
%solve the problem in each tent

for i=1:I_max
    %selection of valid value since -1 is a defualt value 
    fprintf('Solving the problem inside the tents %d ...\n',i);



    time = tents(:,2,i);
    x_ = tents(:,1,i);
    acceptable_x = x_(x_>-1);
    acceptable_time = time(time>=0);
    delta = max(acceptable_time)-min(acceptable_time);
    %mapp the physical tent Ki into the simplest tent Ki_hat
    [front,L_] = compute_analytical_front1(acceptable_x,Dati,acceptable_time);
    L = L_;
    %solve the problem in Ki_hat 
    [U_hat_i,V_hat_i] = solver(Dati,tentDependencies(i,:),U,V,acceptable_x,acceptable_time,front);
    
    if Dati.visual_graph_3D == 'Y'
        surface_tent_solution(acceptable_x,front,U_hat_i,V_hat_i,Dati)
    end

    U{end+1} = U_hat_i;
    V{end+1} = V_hat_i;
    %t = 0:Dati.dt:1;
    
    fprintf('mapping back the solution of the tent %d...\n',i);

    [u,v] = map_solution(U_hat_i,V_hat_i,front,acceptable_x,acceptable_time,Dati);
    sol_u{end+1} = u;
    sol_v{end+1} = v;
    if i>= (I_max-n_even +1)
        if length(acceptable_x) ==4
            list_u = [list_u; U_hat_i(2:3,end)];
            list_v = [list_v;V_hat_i(2:3,end)];
        else
            if sum(abs(acceptable_x-x0)<1e-5)>0
                T = max(acceptable_time);
                list_u = [list_u;U_hat_i(:,end)];
                list_v = [list_v;V_hat_i(:,end)];
            else
                list_u = [list_u;U_hat_i(2,end)];
                list_v = [list_v;V_hat_i(2,end)];
            end
        
        end
    end
    
    
    
end

if Dati.visual_graph == 'Y'
    general_plot(n_even,n_odd,sol_u,sol_v,Dati,nodes,delta)
end

fprintf('evaluating the error...\n');


error = compute_error1(Dati,list_u,list_v,T);

end