function [cs3avg,cs1avg] = get_avgConcentration(cs3,cs1,data,matrices_spm,nodes)
% GET_AVGCONCENTRATION This function calculates the average lithium 
% concentration in each particle.
% 
% INPUTS
% cs3           Concentration at cathode [mol/m3], see get_concentration
% cs1           Concentration at anode [mol/m3], see get_concentration
% data          Structure containg the model parameters, see get_modelData
% matrices_spm  Structure containing the model matrices, see get_model
% nodes         Structure containing the Chebyshev nodes, see get_nodes
%
% OUTPUTS
% cs3avg        Average concentration in cathode [mol/m3]
% cs1avg        Average concentration in anode [mol/m3]
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

xr1 = nodes.xc2xp(nodes.xr,'r1')*ones(1,size(cs1,2));
xr3 = nodes.xc2xp(nodes.xr,'r3')*ones(1,size(cs3,2));

cs1avg = ((3/data.Rs1^3)*data.Rs1*matrices_spm.wn)*(xr1.^2.*cs1)/2;
cs3avg = ((3/data.Rs3^3)*data.Rs3*matrices_spm.wn)*(xr3.^2.*cs3)/2;