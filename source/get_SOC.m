function [SOC1,SOC3] = get_SOC(cs1_avg,cs3_avg,data)
%GET_SOC This function computes the bulk SOC of the anode SOC1 and the
%cathode SOC3.
%
% INPUTS
% cs1_avg       Average anode concentration [mol/m3]
% cs3_avg       Average cathode concentration [mol/m3]
% data          Structure containg the model parameters, see get_modelData
% 
% OUTPUTS
% SOC3          Bulk SOC based on the cathode average concentration.
% SOC1         	Bulk SOC based on the anode average concentration.
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

SOC3 = 100*(cs3_avg/data.cs3_max - data.y3_soc0)/...
                        (data.y3_soc1 - data.y3_soc0);
SOC1 = 100*(cs1_avg/data.cs1_max - data.x1_soc0)/...
                        (data.x1_soc1 - data.x1_soc0);