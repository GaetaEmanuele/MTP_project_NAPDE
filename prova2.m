clear all
close all
clc

Dati = C_dati();
h = Dati.h;
h_vec = [h,h/2,h/4];
Dati.visual_graph = 'N';
Dati.visual_graph_3D = 'N';

error_vec = zeros(size(h_vec));
for i = 1:length(h_vec)
    Dati.h = h_vec(i);
    [~,~,~,~,err] = c_main(Dati);
    error_vec(i) = err;
end

figure
loglog(h_vec,error_vec)
hold on 
loglog(h_vec,h_vec)
loglog(h_vec,h_vec.^2)
legend('error','h','h^2')
