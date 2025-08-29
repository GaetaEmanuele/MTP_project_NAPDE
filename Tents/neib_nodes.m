function Nl = neib_nodes(x)
    N = length(x);
    Nl = zeros(2,N);
    Nl(1,1)=-1;
    Nl(2,1) = 2;
    Nl(1,end) = N-1;
    Nl(2,end) = -1;
    for i=2:(N-1)
        Nl(1,i) = i-1;
        Nl(2,i) = i+1;
    end
end