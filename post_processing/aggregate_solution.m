function [index_i,index_j] = aggregate_solution(nodes,x,t,Dati,N)
    h = Dati.h;
    index_i = zeros(length(t),1);
    index_j = min(floor((x - Dati.domain(1)) / h) + 1, N);
    j = 2;
    bool = false;
    tol = 1e-4;
    for i=1:length(t)
        while ~bool
            for k = 1:size(nodes,2)
                if abs(nodes(j,k)-t(i))<tol
                    index_i(i) = j - 1; %becuse the 1st row of nodes are xi
                    bool = true;
                end
            end
            j = j+1;
        end
        bool = false;
        j = 2;
    end
end