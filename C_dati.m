%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 
%                'Domain',            % set the domain [x1,x2]
%                'rho',               % set the parameter rho
%                'mu',                % set the parameter mu
%                'exact_sol',         % set the exact solution 
%                'force',             % set the forcing term in space
%                'force_time',        % set the forcing term in time
%                'derivativet',       % v exact
%                'derivativex',       % u exact
%                'h',                 % discrtitization parameter in space
%                'dt',                % local dt for each tent
%                'T',                 % max time considered
%                'element',           % elment used, defualt on P1 
%                'visual_graph',      % if you want to display the graphical results ['Y','N']
%                'visual_graph_3D',   % if you want to have 3D graph of the solution for each tent ['Y','N']
%                'plot_errors'        % you want to plot the computed errors ['Y','N']
% 
%========================================================================================================

function [Dati]=C_dati()
Dati = struct(  'domain',           [0,1],... 
               ... % Domain bounds
               'rho',               1,...
               ... % rho
               'mu',               1,...
               ... % mu
               'exact_sol',        'sin(pi*x).*cos(pi*t)',...      
               ... % Definition of exact solution
               'force',            '0.*x',...  
               ... % Forcing term
               'force_time',     '0.*t',...    
               ... % Definition of exact derivative 
               'derivativet',     '-pi*sin(pi*x).*sin(pi*t)',...    
               ... % Definition of exact derivative
               'derivativex',     'pi *cos(pi*x).*cos(pi*t)',...    
               ... % Definition of exact derivative
               'h',                 0.05,...   
               ... % step in space
               'dt',                  0.005,...   
               ... % step in space
               'T',                 0.3,...   
               ... % Max time
               'element',                 'P1',...   
               ... % element P1 or P2
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'visual_graph_3D',      'N',...
               ... % Visualization of the solution inside the single tent
               'plot_errors',       'Y' ...
               ...% Compute Errors
               );
end



