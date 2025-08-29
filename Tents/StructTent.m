function Tent = StructTent(x_coordinate, t_higher)
    % Definition of the data structure of a single tent
    A = zeros(4,2);
    if length(t_higher)==3
        A(1:3,1) = x_coordinate;
        A(1:3,2) = t_higher;
        A(4,1) = -1;
        A(4,2) = -1;
    elseif length(t_higher)==4
        A(:,1) = x_coordinate;
        A(:,2) = t_higher;
    else
        disp('error dimension of dimensions')
    end
    Tent = A;
end

