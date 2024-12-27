function Gparray = Gppfun(l,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lstar = param.lstar;
a = param.a;
Gmult = param.Gmult;

lscaled = l/lstar;

Gparray = Gmult*(2*a/lstar^2-a*lscaled/lstar^2+1/lstar^2).*exp(-lscaled);

end