function [ postproc_result ] = get_postproc( result,data,nodes,matrices_spm,I)
%GET_POSTPROC This function postprocess the time/state results obtained
%from the ODE solver to compute the voltage, temperature, current,
%concentration cs, average concentration and anode/cathode
%concentration-based SOC.
%
% INPUTS
% result        Structure containing a vector of time and corresponding
%               model states returned by the ODE solver.
% data          Structure containg the model parameters, see get_modelData
% nodes         Structure containing the Chebyshev nodes, see get_nodes
% matrices_spm  Structure containing the model matrices, see get_model
% I             Handle to an anonymous function returning the input current
%               as a function of time.
%
% OUPUTS
% postproc_result   A structure containing the 'original' result structure
%                   and the computed values of voltage, temperature,
%                   concentration profiles, SOC,...
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

postproc_result.time    = real(result.time);        % Time [s]
postproc_result.state   = real(result.state);    	% Model states
%{
    We are taking the real part only of the time and state as sometimes the
    ODE solver finds complex solutions with infinitely small imaginary
    parts. This would mess up the plotting of the results.
%}  

postproc_result.current = I(postproc_result.time);              % Input current [A]
postproc_result.i_app   = postproc_result.current/data.As;      % Current density [A.m-2]
postproc_result.j3      = ...
    -postproc_result.i_app/data.as3/data.F/data.thick3;         % Cathode reaction rate [mol.m-2.s-1]
postproc_result.j1      = ...
     postproc_result.i_app/data.as1/data.F/data.thick1;         % Anode reaction rate [mol.m-2.s-1]

% Computing the voltage and temperature
[V,T] = get_measurements(postproc_result.time',...
    postproc_result.state',I,data,matrices_spm);
postproc_result.voltage = V';
postproc_result.temperature = T';

% Computing the concentration profile in anode cs1 and cathode cs3 (value 
% at the Chebyshev discretization nodes.
[cs3,cs1] = get_concentration(postproc_result.time',postproc_result.state',I,data,matrices_spm);
postproc_result.cs1 = cs1';
postproc_result.cs3 = cs3';

% Computing the average concentration in anode cs3avg and cathode cs1avg
[cs3avg,cs1avg] = get_avgConcentration(cs3,cs1,data,matrices_spm,nodes);
postproc_result.cs1avg = cs1avg';
postproc_result.cs3avg = cs3avg';

% Computing the concentration-based SOC in anode SOC1 and cathode SOC3
[SOC1,SOC3] = get_SOC(cs1avg,cs3avg,data);
postproc_result.SOC1 = SOC1';
postproc_result.SOC3 = SOC3';