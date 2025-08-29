%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test 
%                'Domain',            % set the domain [x1,x2;y1,y2]
%                'exact_sol',         % set the exact solution
%                'force',             % set the forcing term
%                'grad_exact_1',      % set the first componenet of the gradient of the exact solution
%                'grad_exact_2',      % set the second componenet of the gradient of the exact solution
%                'fem',               % set finite element space
%                'nqn_1D',            % number of quadrature nodes for integrals over lines
%                'nqn_2D',            % number of quadrature nodes for integrals over elements
%                'MeshType',          % set the elements of the mesh  'TS'
%                'refinement_vector', % set the level of refinement for the grid
%                'visual_graph',      % if you want to display the graphical results ['Y','N']
%                'print_out',         % if you want to print out the results ['Y','N']
%                'plot_errors'        % you want to plot the computed errors ['Y','N']
% 
%========================================================================================================

function [Dati]=C_dati()
Dati = struct(  'domain',           [0,1],... 
               ... % Domain bounds
               'c',               1,...
               ... % Velocity of propagation
               'exact_sol',        'sin(pi*x)',...      
               ... % Definition of exact solution
               'force',            'pi*pi*sin(pi*x)',...  
               ... % Forcing term
               'force_time',     '0.*t',...    
               ... % Definition of exact derivative 
               'derivative',     'pi*cos(pi*x).*cos(pi*t)',...    
               ... % Definition of exact derivative
               'boundary_cond',     '0.*x',...    
               ... % boundary condition space
               'boundary_cond_time',     '0.*t',...    
               ... % boundary condition time
               'initial_cond',     'sin(pi*x)',...    
               ... % initial condition
               'initial_velocity',     '0.*x',...    
               ... % initial condition
               'h',                 0.1,...   
               ... % step in space
               'dt',                 0.005,...   
               ... % step in space
               'T',                 2,...   
               ... % Max time
               'element',                 'P1',...   
               ... % element P1 or P2
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors
               );
end



