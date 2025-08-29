function [bd_cond] = check_bd_condition(x_,x0,xN)
    %the function create a struct of 2 fields 
    %where we save the type of bd condition for the ending nodes
    %Dirichlet for physical nodes Sommerfeld for artifial boundary
    if x_(1) == x0
        bd_cond.field1 = 'Dirichlet';
    else 
        bd_cond.field1 = 'None'
    end

    if x_(end) == xN
        bd_cond.field2 = 'Dirichlet';
    else 
        bd_cond.field2 = 'None'
    end

end