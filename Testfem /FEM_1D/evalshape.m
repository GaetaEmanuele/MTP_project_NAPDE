function [dphiq,Deriv] = evalshape(basis,node)
%% [dphiq,Grad] = C_evalshape(basis,node_2D)
%==========================================================================
% Evaluation of basis functions and gradients on quadrature nodes 
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          basis    : (struct)  see C_shape_basis.m
%          node2D   : (array real) [nqn x 2] (x,y) coord of quadrature nodes   
%
%    OUTPUT:
%          dphi     : dphiq(:,:,i) vecotr containing the evaluation of the
%                     i-th basis function on the quadrature nodes 
%          Grad     : dphiq(:,1:2,i) vecotr containing the evaluation of the
%                     gradinet of dphi(:,:,i) on the quadrature nodes 


nln = length(basis);

%==========================================================================
% EVALUATION OF BASIS FUNCTIONS AND DERIVATIVES
%==========================================================================

for s=1:nln
    % coordinates
    csi = node;
    % evaluation of basis functions
    dphiq(:,s) = eval(basis(s).fbases);
    % evaluation of gradients
    Deriv(:,s) = eval(basis(s).derivative);
end
