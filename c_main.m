function [sol_u,sol_v,nodes,n_even,error] = c_main(Dati)

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
%creation of the space-time meshes
[tents,nodes,tentDependencies,n_even,n_odd] = create_tents1(Dati);

%N = (Dati.domain(2)-Dati.domain(1))/Dati.h;
%x = linspace(Dati.domain(1),Dati.domain(2),N);
U = {};
V = {};
if Dati.visual_graph == 'Y'
    C_plot(nodes,Dati);
end
I_max = size(tents,3);
n_sol = size(nodes,1)-3;
n_col = size(nodes,2);
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
    
    U{end+1} = U_hat_i;
    V{end+1} = V_hat_i;
    %t = 0:Dati.dt:1;
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

error = compute_error1(Dati,list_u,list_v,T);

end