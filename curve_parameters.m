function param = curve_parameters
%
% Called by the program curve_evolution
% sets various parameters for that program, including the spline 
% points of the initial curve
%
% Set 1 to make a movie
param.movie = 0;
param.movie_name = "curve1_hires.avi";

% Set to 1 if you want the area to be preserved with a normal 
%    velocity projection;
param.conserving = 0;
%
% Gcurve parameters
param.lstar = 0.25;
param.a = 1;
param.Gmult = 10; % Set to zero for no adhesion terms 
%
% Normal velocity f parameters
param.epsilon=0.1;
param.mu=1;
param.rho=1.5;

% Maximum number of Newton iterations
param.maxNewton=100;
% Tolerance of Newton iterations
param.tol=1e-8;
% Maximum domain size for Newton iterates of the curve positions
param.maxSize = 10;
% Error tolerance per time step
param.sigma = 1e-4;

% First time step 
param.firstk = 1e-5;
% Maximum fractional time step increase
param.alpha = 1.3;
% Minimum time step for failure
param.mink = 1e-12;
% Maximum time step
param.maxk = 1;

return
