function tentDependencies = findTentDependencies1(tents,n_even,n_odd)
    %inputs n_even and n_nodes says how much tents are built in
    %even and odd iteration 

    
    numTents = size(tents,3);  % total number of tent

    % Creation of a matrix of dependence
    % all the row associated with the jth column
    % are the tents that tent j must wait to be used
    tentDependencies = zeros(numTents,2);
    i=0;
    ii=1;
    while ii<numTents
        if (mod(i,2)==0) && (i>0)
            tentDependencies(ii,1) = ii - n_odd;
            for j=(ii+1):(ii+n_even-2)
                tentDependencies(j,1) = j-n_odd-1;
                tentDependencies(j,2) = j-n_odd;
            end
            tentDependencies(ii+n_even-1,1) = ii+n_even-2 - n_odd;
            ii = ii + n_even;
        elseif mod(i,2)>0
            for j=ii:(ii+n_odd-1)
                tentDependencies(j,1) = j-n_even;
                tentDependencies(j,2) = j-n_even+1;
            end
            ii = ii + n_odd;
        elseif i ==0
            ii = ii + n_even;
        end
        i = i+1;
    end

end