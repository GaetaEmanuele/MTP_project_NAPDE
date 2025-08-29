function coord = map(connectivity,xq,x)
x0 = x(connectivity(1));
x1 = x(connectivity(2));
coord = (x1 - x0) / 2 * xq + (x0 + x1) / 2 ;
end