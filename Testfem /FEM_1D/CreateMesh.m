function [ L, lnn, nln, el, connectivity, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  )
% Creates mesh of Line with 1D linear and quadratic elements. 
% Node numbering
%     P1:         1*-------------*2
%     P2:         1*------2------*3
% Important outputs:
%     egnn:       Element Global Node Numbering       (tne x lnn)  Connectivity Matrix
%     x:          Nodal Coordinates Vector
% 
% Length of the domain: L
L = xend - xstart;

if strcmp(elementtype, 'P1')
    % Local node numbering
    lnn = 1 : 2;
    % Number of nodes in the element
    nln = size(lnn,2);
    % Length of the element.
    el = L / tne;
    % Total number of nodes in the mesh
    tnn = (tne*(size(lnn,2)-1))+1;
    connectivity = zeros(size(lnn,2),tne);
    % Connectivity Matrix (egnn = Element Global Node Numbering)
    i = 1;
    for en = 1 : tne
        for I = lnn
            connectivity(I,en) = i;
            i = i + 1;
        end
        i = i - 1;
    end
    % Nodal Coordinates
    x  = (xstart : el : xend)';
    
elseif strcmp(elementtype, 'P2')
    % Local node numbering
    lnn = 1 : 3;
    % Number of nodes in the element
    nln = size(lnn,2);
    % Length of the element.
    el = L / tne;
    % Total number of nodes in the mesh
    tnn = (tne*(size(lnn,2)-1))+1;
    connectivity = zeros(size(lnn,2),tne);
    % Connectivity Matrix (egnn = Element Global Node Numbering)
    i = 1;
    for en = 1 : tne
        for I = lnn
            connectivity(I, en) = i;
            i = i + 1;
        end
        i = i - 1;
    end
    % Nodal Coordinates
    x  = (xstart : el/2 : xend)';
    
end