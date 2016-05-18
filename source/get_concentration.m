function [cs3, cs1] = get_concentration(t,y,I,data,matrices_spm)
% GET_CONCENTRATION Function to calculate the concentration at the 
% different Chebyshev nodes in each particle.
%
% INPUTS
% t             Time in seconds at which the concentration has to be 
%               calculated.
% y             The state at which the concentration has to be calculated.
% I             Handle to an anonymous function returning the input current
%               as a function of time.
% data          Structure containg the model parameters, see get_modelData.
% matrices_spm  Structure containing the model matrices, see get_model.
%
% OUTPUTS
% cs3           Concentration at cathode [mol/m3]. 
%                   cs3(N+1) = concentration at node at centre of sphere (r = 0)
%                   cs3(1)   = concentration at node on surface of sphere
% cs3           Concentration at anode [mol/m3]. 
%                   cs3(N+1) = concentration at node at centre of sphere (r = 0)
%                   cs3(1)   = concentration at node on surface of sphere
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

    N = matrices_spm.N;
    M = matrices_spm.M;

    % extract parameters from the state
    us3 = y(1:N-1,:);           % u at cathode
    us1 = y(N:2*N-2,:);         % u at anode
    T   = y(2*N-1,:);           % battery temperature
  
    i_app = I(t)/data.As;   % Current density [A.m-2]
    
    % calculate molar flux on the electrodes [mol.m-2.s-1]
    j3 = -i_app/data.as3/data.F/data.thick3;   	% on cathode
    j1 =  i_app/data.as1/data.F/data.thick1;   	% on anode
    
    % Calculate the temperature-dependent diffusion coefficients according
    % to an Arrhenius relation D = Dref * exp (E/R * (1/Tref - 1/T))
        % Dref = D at reference temperature
        % E = activation energy
        % R = gas constant
        % Tref = temperature at which the reference D was measured
        % T = temperature at which we want to calculate D
    Ds3_var = data.Ds3_ref*exp(data.Ea_Ds3/data.R*(1/data.T_ref - 1./T));   % D at cathode
    Ds1_var = data.Ds1_ref*exp(data.Ea_Ds1/data.R*(1/data.T_ref - 1./T));   % D at anode
    
    % Calculate connecntration of inner nodes (=all nodes except the
    % centre node at r=0) using predefined matrix
    cs3_inner  = matrices_spm.C3*us3 + ...  % cathode concentration
        (matrices_spm.D3*j3)./(ones(N,1)*Ds3_var);
    cs1_inner  = matrices_spm.C1*us1 + ...  % anode concentration
        (matrices_spm.D1*j1)./(ones(N,1)*Ds1_var);
    
    % calculate concentration at the centre node (r=0) using the boundary
    % condition at the particle surface
    % @x=1 dc/dx = -j*Rs^2/D_s
    
    cs3_center = (-1/matrices_spm.DM(1,N+1,1))*(...
        (matrices_spm.DM(1,1:N,1) + matrices_spm.DM(1,N+2:M+1,1)*matrices_spm.P )*cs3_inner + ...
        j3*data.Rs3./Ds3_var );
    cs1_center = (-1/matrices_spm.DM(1,N+1,1))*(...
        (matrices_spm.DM(1,1:N,1) + matrices_spm.DM(1,N+2:M+1,1)*matrices_spm.P )*cs1_inner + ...
        j1*data.Rs1./Ds1_var );

    % Combining the results [c_inner;c_center]
    cs3 = vertcat(cs3_inner,cs3_center);
    cs1 = vertcat(cs1_inner,cs1_center);
end
    