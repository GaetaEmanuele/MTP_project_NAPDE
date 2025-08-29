function [error] = compute_error(Dati,u,v)


u = cell2mat(u);
v = cell2mat(v);

u_ex = inline(Dati.derivativex,'x','t');
v_ex = inline(Dati.derivativex,'x','t');

x0 = Dati.domain(1);
x1 = Dati.domain(2);
t = Dati.T;
h = Dati.h;
x_temporary = x0:2*h:x1;
tne = length(x_temporary)-1;

[ ~, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( Dati.element, tne, x0, x1);
    
[basis] = shape_basis(Dati.element);

[xq, W_1D] = Quadrature(nln);
xq = sort(xq);

[dphiq,Deriv] = evalshape(basis,xq);
M = sparse(tnn,tnn);
nq = length(W_1D);

for ie=1:tne
        iglo = connectivity(1:nln,ie);
        M_loc = zeros(nln,nln);
        
        for q=1:nq
            %B = BJ(q);
            for i=1:nln
                deriv_i = Deriv(1,i);
                phi_i = dphiq(i,q);
                
                for j=1:nln
                    deriv_j = Deriv(1,j);
                    phi_j = dphiq(j,q);
                    
                    M_loc(i,j) = M_loc(i,j) + (h)*(phi_i*phi_j)*W_1D(q);
                    
                end
            end
        end
        M(iglo,iglo) = M(iglo,iglo) + M_loc;
        
        
end

uexa = u_ex(x,t);
vexa = v_ex(x,t);
diff1 = uexa -u';
diff2 = vexa -v';
err1 = (diff1)' * M * diff1;
err2 = (diff2)' * M * diff2;
error = sqrt(err1+err2);


end


