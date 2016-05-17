function [ initSPM ] = get_init( x1_init,y3_init,T_init,data,nodes,model )
%GET_INIT This function compute the initial state vector based on the
% user-defined initial conditions (initial anode and cathode stoichiometry,
% assumed uniform within each particle, amd initial temperature).
%
% INPUTS
% x1_init       Initial anode stoichiometry
% y3_init       Initial cathode stoichiometry
% T_init        Initial battery temperature
% data          Structure containg the model parameters, see get_modelData
% nodes         Structure containing the Chebyshev nodes, see get_nodes
% model         Structure containing the model matrices, see get_model
%
% OUTPUTS
% initSPM       Structure containing the initial state of the battery, in
%               particular the initial state vector initSPM.y0  
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

N = model.N;

initSPM.cs1 = data.cs1_max*x1_init;   % Initial anode concentration
initSPM.cs3 = data.cs3_max*y3_init;   % Initial cathode concentration
initSPM.us1 = initSPM.cs1.*nodes.xc2xp(nodes.xr(2:N),'r1');
initSPM.us3 = initSPM.cs3.*nodes.xc2xp(nodes.xr(2:N),'r3');
initSPM.T = T_init;

% Initial state vector to be provided to the ODE solver
initSPM.y0 = [initSPM.us3;initSPM.us1;initSPM.T];

end

