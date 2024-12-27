% Script to generate the Figure in the Numerical Approximation section of
% the paper, "GRADIENT FLOWS OF INTERFACES: ADHESION AND INEXTENSIBLE
% ADVECTIVE TRANSPORT" by Promislow and Wetton.

% Generate the initial data object state0
% Edit pcurve_init_v2 to change the number of spatial grid points N
fprintf('*** Generate the initial data in data object state0 ***\n \n') 
pcurve_init_v2;
state0

fprintf('*** Hit any ket to continue ***\n \n') 
pause

% Compute to time 0.1 with the result in state1
% Edit curve_parameters to change the local error tolerance sigma
fprintf('*** Compute to time 0.1 with the result in state1 ***\n \n') 

state1 = curve(state0, 0.1);
state1


