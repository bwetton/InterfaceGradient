function [x,y,dxs,dys] = cinit_xy(sigma)

% Ellipse 
% ellipticity = 4;
% theta = 2*pi*sigma;
% eps=0.10;
% 
% x = cos(theta)+eps*(sin(5*theta)).^2;
% y = (sin(theta)+eps*cos(3*theta))/ellipticity;
% dxs = 2*pi*(-sin(theta)+2*eps*5*cos(5*theta).*sin(5*theta));
% dys = 2*pi*(cos(theta)-3*eps*sin(3*theta))/ellipticity;

% Base Trillium curve
theta = 2*pi*sigma;
eps = 0.35;

r = 1 + eps*cos(3*theta);
dxr = - eps*3*sin(3*theta)*2*pi;
x = r*cos(theta);
y = r*sin(theta);
dxs = dxr*cos(theta) - r*2*pi*sin(theta);
dys = dxr*sin(theta) + r*2*pi*cos(theta);

% Add a perturbation 
x = x + 0.05*sin(4*theta);
dxs = dxs + 0.05*4*2*pi*cos(4*theta);

end

