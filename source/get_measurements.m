function [V,T] = get_measurements(t,y,I,data,matrices_spm)
%GET_MEASUREMENTS This function computes the battery terminal voltage and
% temperature from the state vector y.
%
% INPUTS
% t             Time in seconds at which to compute voltage and temperature
% y             State of the model
% I             Handle to an anonymous function returning the input current
%               as a function of time.
% data          Structure containg the model parameters, see get_modelData
% matrices_spm  Structure containing the model matrices, see get_model
%
% OUTPUTS
% V             Battery terminal voltage in volts.
% T             Battery temperature in Kelvins
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

% For comments, see function derivs_spm.m

    N = matrices_spm.N;
    
    i_app = I(t)/data.As;       % Applied current density
    
    us3 = y(1:N-1,:);
    us1 = y(N:2*N-2,:);
    
    T   = y(2*N-1,:);
    
    j3 = -i_app/data.as3/data.F/data.thick3;
    j1 =  i_app/data.as1/data.F/data.thick1;

    Ds3_var = data.Ds3_ref*exp(data.Ea_Ds3/data.R*(1/data.T_ref - 1./T));
    Ds1_var = data.Ds1_ref*exp(data.Ea_Ds1/data.R*(1/data.T_ref - 1./T));
    
    cs3 = matrices_spm.C3*us3 + (matrices_spm.D3*j3)./(ones(N,1)*Ds3_var);
    cs1 = matrices_spm.C1*us1 + (matrices_spm.D1*j1)./(ones(N,1)*Ds1_var);

%     css3 = cs3(N,:);
%     css1 = cs1(N,:);
    css3 = cs3(1,:);
    css1 = cs1(1,:);

        k3_T    = data.k3_ref*exp(data.Ea_k3/data.R*(1/data.T_ref - 1./T));
        k1_T    = data.k1_ref*exp(data.Ea_k1/data.R*(1/data.T_ref - 1./T));
        i0_3    = k3_T*data.F.*(data.ce_avg.^0.5).*(css3.^0.5).*(data.cs3_max-css3).^0.5;
        i0_1    = k1_T*data.F.*(data.ce_avg.^0.5).*(css1.^0.5).*(data.cs1_max-css1).^0.5;
    
    eta3 = (2*data.R*T/data.F).*asinh(-0.5*i_app/data.as3/data.thick3./i0_3);
    eta1 = (2*data.R*T/data.F).*asinh( 0.5*i_app/data.as1/data.thick1./i0_1);
        
        % Particle surface stoichiometry
        x3 = css3/data.cs3_max;
        x1 = css1/data.cs1_max;
        [ U_ocp1,~,U_ocp3,~] = get_openCircuitPotential( x1,x3,T,data );
        
        % Cell terminal voltage
        V   =    (U_ocp3 + eta3 ) - (U_ocp1 + eta1 ) - data.Rc*i_app;        
end