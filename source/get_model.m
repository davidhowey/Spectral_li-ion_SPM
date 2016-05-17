function model = get_model( data,nodes,N )
%GET_MODEL This function computes the matrices A, B, C and D of the 
% particel diffusion state space model of the single particle model for
% each electrode:
%    du/dt  = A*u + B*j
%    cs     = C*u + D*j
% with:
%   u=rc    state of the diffusion model.
%   j       volumetric reaction rate.
%   cs      particle concentration profile (including at the surface)
%
% INPUTS
% data          Structure containg the model parameters, see get_modelData
% nodes         Structure containing the Chebyshev nodes, see get_nodes
% N             The number N+1 of Chebyshev nodes used to discretize the
%               diffusion PDE in each particle.
%
% OUTPUTS
% model         A structure constaining the matrices A, B, C and D for each
%               electrode, and also the Chebyshev differentiation matrices
%               DNr and the Clenshaw-Curtis quadrature weights wn.
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.


% M+1: number of Chebyshev nodes in each particle on the domain [-R,R]
% mapped onto [-1,1], with N = M/2 defined by the user.
M = 2*N;
model.N = N;
model.M = M;

% Computing the Chebyshev differentiation matrices
[~,DM] = chebdif(M+1,2);
model.DM = DM;
DM1 = squeeze(DM(1,:,1));
DM2 = squeeze(DM(1:N,:,2));

P   = fliplr(eye(N)); model.P = P;  % Permutation matrix, backward-identity

% Modified differentation matrices accounting for the solution symmetry
DN2 = DM2(:,1:N) - DM2(:,N+2:M+1)*P;
DN1 = DM1(1,1:N) - DM1(1,N+2:M+1)*P;


A = DN2(2:N,2:N) + (DN2(2:N,1)*DN1(1,2:N))/(1-DN1(1,1));
B = DN2(2:N,1)/(1-DN1(1,1));
C = DN1(1,2:N)/(1-DN1(1,1));
D = 1/(1-DN1(1,1));

% Anode diffusion state-space model
model.A1 = (1/data.Rs1^2)*A;
model.B1 = B;
model.C1 = vertcat( ...
            C/data.Rs1 , ...
            diag(1./nodes.xc2xp(nodes.xm(2:N),'r1')) );
model.D1 = vertcat(...
            data.Rs1*D , ...
            zeros(N-1,1) );

% Cathode diffusion state-space model
model.A3 = (1/data.Rs3^2)*A;
model.B3 = B;
model.C3 = vertcat( ...
            C/data.Rs3 , ...
            diag(1./nodes.xc2xp(nodes.xr(2:N),'r3')) );
model.D3 = vertcat( ...
           	data.Rs3*D , ...
            zeros(N-1,1) );

% Clenshaw-Curtis quadrature weights for calculating integrals
[~,wm] = clencurt(M);
wn = horzcat( ...
    (wm(1:N) + wm(N+2:M+1)*P) , ...
    wm(N+1) );
model.wn = wn;

end

