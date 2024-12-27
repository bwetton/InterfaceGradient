function [f, dL, dk, dkss] = curvef(L, L0, kappa, kappass,param)

% L0 has the initial length of the curve if needed
% f is the normal velocity, [dL dk dkss] its partial derivatives

mu = param.mu;
rho = param.rho;
epsilon = param.epsilon;

f = -epsilon*kappass - epsilon*kappa.^3/2 +mu*(L-rho*L0)*kappa; 

dL = kappa;
dk = -3/2*epsilon*kappa.^2+mu*(L-rho*L0);
dkss = -epsilon;

