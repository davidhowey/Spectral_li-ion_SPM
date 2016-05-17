function [ V1,dV1dT,V3,dV3dT] = get_openCircuitPotential( x1,x3,T,data )
% GET_OPENCIRCUITPOTENTIAL This function computes the open-circuit
% potential and entropy coefficient of the anode and cathode as a function
% of lithium stoichiometry and battery temperature.
%
% INPUTS:
% x1            stoichiometry of anode active material.
% x3            stoichiometry of cathode active material.
% T             battery temperature.
%
% OUTPUTS:
% V1            anode open-circuit potential.
% dV1dT         anode entropy coefficient.
% V3            cathode open-circuit potential.
% dV3dT         cathode entropy coefficient.
%
%
%   The LCO OCP and entropy coefficient fitting function provided in this 
%   file are taken from:
%       -   "Lithium-ion cell modeling using orthogonal  collocation on
%       finite elements", Cai & White; JPS 217 (2012) 248-225.
%       -   "Development of First-Principles Capacity fade model for 
%       li-ion cells"; Ramadass et al. JES 151(2) A196-A203 (2004).
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.


    %% ANODE (Graphite - LiC6)
    % Open-circuit potential at reference temperature 25dC
    V1_ref   = ...
        0.7222 + 0.1387*x1 + 0.0290*x1.^(1/2) - 0.0172./x1 + ...
        0.0019./(x1.^(1.5)) + 0.2808*exp(0.90-15*x1) - ...
        0.7984*exp(0.4465*x1-0.4108);
    
    % Entropy coefficient
    c1_num = [-16515.05308;38379.18127;-37147.89470;19329.75490;...
        -5812.27813;1004.91101;-91.79326;3.29927;0.00527];
    c1_den = [165705.85970;-385821.16070;374577.31520;-195881.64880;...
        59431.30001;-10481.80419;1017.23480;-48.09287;1];
    dV1_dT_num = c1_num(1);
    dV1_dT_den = c1_den(1);
    for i = 2:7
        dV1_dT_num = dV1_dT_num.*x1 + c1_num(i);
        dV1_dT_den = dV1_dT_den.*x1 + c1_den(i);
    end
    dV1dT = 1e-3*dV1_dT_num./dV1_dT_den;
    
    % Open-circuit potential at temperature T
    V1 = V1_ref + (ones(size(x1,1),1)*T-data.T_ref).*dV1dT;

    %% CATHODE (Cobalt oxide - LiCoO2)
    % Open-circuit potential at reference temperature 25dC
    V3_ref   = ...
        ( -4.656 + 88.669*x3.^2 - 401.119*x3.^4 + 342.909*x3.^6 - ...
           462.471*x3.^8 + 433.434*x3.^10)./...
        ( -1 + 18.933*x3.^2 - 79.532*x3.^4 + 37.311*x3.^6 - ...
           73.083*x3.^8 + 95.96*x3.^10);
    
    % Entropy coefficient
    c3_num = [0;0.61154;-1.36455;0.92837;-0.19952];
    c3_den = [3.04876;-9.82431;11.47636;-5.66148;1];
    dV3dT_num = c3_num(1);
    dV3dT_den = c3_den(1);
    for i = 2:5
        dV3dT_num = dV3dT_num.*x3 + c3_num(i);
        dV3dT_den = dV3dT_den.*x3 + c3_den(i);
    end
    dV3dT = 1e-3*dV3dT_num./dV3dT_den;
    
    % Open-circuit potential at temperature T
    V3 = V3_ref + (ones(size(x3,1),1)*T-data.T_ref).*dV3dT;

end

