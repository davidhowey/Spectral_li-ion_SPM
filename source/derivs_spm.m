function dy = derivs_spm(t,y,I,data,matrices_spm)  
% DERIVS_SPM This function computes the time derivative f of the dynamic
% equation of the SPM thermal-electrocehmical battery state-space model
% dy/dt = f(t,y)
%
% INPUTS
% t             Time in seconds at which we want to calculate the time 
%               derivative.
% y             State at which we want to calculate the time derivative.
% I             The applied current [A] (model input) This must be the 
%               handle of an anonymous function returning input current
%               as a function of time t.
% data          Structure with parameters of the battery 
%               (as built by 'model_parameters/get_modelData')
% matrices_spm  Structure with differentiation matrices
%               (as built by 'source/get_model').
%
% OUTPUTS
% dy            Time derivative of the state
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.
    
    N = matrices_spm.N;    
    %% Electrochemical model
    % Calculate the diffusion of Li in the electrodes

    % Extract state variables
    us3 = y(1:N-1,:);           % transformed concentration at cathode
    us1 = y(N:2*N-2,:);         % transformed concentration at anode
    T   = y(2*N-1,:);           % battery temperature
    
    i_app = I(t)/data.As;       % Applied current density [A.m-2]
    
    % calculate molar flux on the electrodes [mol.m-2.s-1]
    j3 = -i_app/data.as3/data.F/data.thick3;         % at cathode
    j1 =  i_app/data.as1/data.F/data.thick1;         % at anode

    % Calculate the temperature-dependent diffusion coefficients according
    % to an Arrhenius relation D = Dref * exp (E/R * (1/Tref - 1/T))
        % Dref = D at reference temperature
        % E = activation energy
        % R = gas constant
        % Tref = temperature at which the reference D was measured
        % T = temperature at which we want to calculate D
    Ds3_var = data.Ds3_ref*exp(data.Ea_Ds3/data.R*(1/data.T_ref - 1./T));   % at cathode
    Ds1_var = data.Ds1_ref*exp(data.Ea_Ds1/data.R*(1/data.T_ref - 1./T));   % at anode

    % Calculate time derivative of transformed concentration using
    % differentiation matrices. The diffusion PDE is a linear equation and
    % hence the time derivative can be expressed as a linear state space
    % equation.
    % du/dt = D*A*u + B*j
        % u = transformed concentration (radius * concentration)
        % D = diffusion coefficient
        % A = differentiation matrix for u
        % B = differentiation matrix for the external variable (current)
        % j = molar flux
    d_us3 = (ones((N-1),1)*Ds3_var).*(matrices_spm.A3*us3) + ... % at cathode
        matrices_spm.B3*j3;
    d_us1 = (ones((N-1),1)*Ds1_var).*(matrices_spm.A1*us1) + ... % at anode
        matrices_spm.B1*j1;

    % Calculate the concentration at the surface of the spheres using
    % predefined matrices
%     css3  = matrices_spm.C3(N,:)*us3 + ...   % at cathode
%         (matrices_spm.D3(N,:)*j3)/Ds3_var;
%     css1  = matrices_spm.C1(N,:)*us1 + ...   % at anode
%         (matrices_spm.D1(N,:)*j1)/Ds1_var;
    css3  = matrices_spm.C3(1,:)*us3 + ...   % at cathode
        (matrices_spm.D3(1,:)*j3)/Ds3_var;
    css1  = matrices_spm.C1(1,:)*us1 + ...   % at anode
        (matrices_spm.D1(1,:)*j1)/Ds1_var;

    % Calculate the temperature-dependent chemical rate constants according
    % to an Arrhenius relation k = kref * exp (E/R * (1/Tref - 1/T))
        % kref = k at reference temperature
        % E = activation energy
        % R = gas constant
        % Tref = temperature at which the reference k was measured
        % T = temperature at which we want to calculate k
    k3_T    = data.k3_ref*exp(data.Ea_k3/data.R*(1/data.T_ref - 1./T)); % at cathode
    k1_T    = data.k1_ref*exp(data.Ea_k1/data.R*(1/data.T_ref - 1./T)); % at anode
    
    % Calculate the exchange current density of the Bulter-Volmer equation
    % i0 = k * F * Ce^alpha_a * Cs_surf^alpha_c * (Cs_max - Cs_surf)^alpha_a
        % k = chemical rate constant of the reaction
        % Ce = concentration in the electrolyte (assumed to be constant in
            % the single particle model because electrolyte is not modelled)
        % alpha_a = transfer coefficient of the anodic reaction 
            % (assumed to be 0.5 for symmetric reaction kinetics)
        % Cs_surf = Li concentration at the surface of the electrode
        % alpha_c = transfer coefficient of the cathodic reaction
            % = 1 - alpha_a = 0.5
        % Cs_max = maximum Li concentration in the electrode
    i0_3    = k3_T*data.F.*(data.ce_avg.^0.5).*(css3.^0.5).*(data.cs3_max-css3).^0.5;  % at cathode
    i0_1    = k1_T*data.F.*(data.ce_avg.^0.5).*(css1.^0.5).*(data.cs1_max-css1).^0.5;  % at anode

    % Calculate the overpotentials at each electrode using the
    % Bulter-Volmer equation. If symmetric reaction kinetics are assumed
    % (alpha_a = alpha_c = 0.5), the BV equation reduces to a hyperbolic
    % sine equation which can be inverted to find the overpotential given
    % the currents. (sinh^(-1) = asinh)
    % eta = 2 * R * T / (n*F) * asinh(x*0.5*i/i0)
        % eta = overpotential [V]
        % R = gas constant
        % T = temperature
        % n = number of electrons involved in reaction = 1
        % F = Faradays constant
        % x = sign of equation (+ or - depending on anodic or cathodic
            % reaction)
        % i = current density on the sphere
            % i     = i_app/thick/as
            % as    = specific surface of the sphere = surface of the sphere
                % (single particle) per unit of electrode volume [m2/m3]
            % i_app = current density on the electrode (current per unit of
                % electrode surface)
            % thick = thickness of the electrode
        % i0 = exchange current density
    eta3 = (2*data.R*T/data.F)*asinh(-0.5*i_app/data.as3/data.thick3/i0_3);   % at cathode
    eta1 = (2*data.R*T/data.F)*asinh( 0.5*i_app/data.as1/data.thick1/i0_1);   % at anode

    %% Thermal model
    % Calculate the change in battery temperature.
    % 3 sources of heat + 1 heat transfer all expressed per unit of 
    % electrode volume:
        % reaction heat
        % reversible heat due to entropy changes
        % ohmic heating
        % heat exchange with the environment
    
    % Reversible heat: q_rev = I_vol*T*dU/dT [W/m3]
        % I_vol = current per unit of battery volume
            % = i/V = i/(A*L) = i_app / L
            % i = total current
            % A = surface
            % L = thickness of the battery
        % T = battery temperature
        % dU/dT = entropic coefficient of the reaction at the current Li
        %   concentration (U = open circuit potential)
            % dU/dT = d(Uc-Ua)/dT = dUc/dT - dUa/dT
            % Uc = open circuit potential of the cathode
            % Ua = open circuit potential of the anode
    
        % Entropy change coefficients in [V/K]
        x3 = css3/data.cs3_max;     % current Li fraction in cathode
        x1 = css1/data.cs1_max;     % current Li fraction in anode       
        [ ~,dOCP1_dT,~,dOCP3_dT] = get_openCircuitPotential( x1,x3,T,data );
        
    % Heat exchange with the environment: assume convective heat transfer
    % q_conv = - h * A * (T - T_amb)
        % q_conv = volumetric convective heat transfer rate [W/m3]
        % h = convective heat transfer coefficient [W/(m2*K)]
        % A = battery surface per unit of battery volume [m2/m3]
        % T = battery temperature [K]
        % T_amb = temperature of the environment [K]
    Q_conv  = - (data.h*data.SA_V*(T-data.T_amb));
    
    % total heat generation in the battery = reversible + reaction + ohmic
    % q_tot = q_rev + q_rea + Qc
    % q_tot = I_vol*T*dU/dT + j*(eta_a - eta_c) + i_app^2*Rc*As/Vc
        % q_tot = volumetric heat generation [W/m3]
        % I_vol = current per unit of battery volume = i_app/L
        %  T = battery temperature
        % dU/dT = entropic coefficient of the reaction at the current Li
            % concentration (U = open circuit potential)
        % eta_a = overpotential at the anode
        % eta_c = overpoential at the cathode
    Q_gen   = (-i_app.*T.*(dOCP3_dT - dOCP1_dT))/data.L + ...
              ( i_app.*(eta1 - eta3))/data.L + ...
                i_app.^2*data.Rc*data.As/data.Vc;

    % PDE for battery temperature
    % rho * cp * dT/dt = q_tot - q_conv
        % rho = density of the battery [kg/m3]
        % cp = heat capacity of the battery [J/(kg*K)]
        % dT/dt = time derivative of temperature
        % q_tot = total volumetric heat generation
        % q_conv = convective heat transfer to the environment
    d_T = (1/data.rho/data.Cp)*(Q_gen + Q_conv);

    %% assemble
    % assemble time derivative of the state
    dy = [ d_us3 ; d_us1 ; d_T ];
    
end