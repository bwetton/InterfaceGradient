function Garray = Gfun(l,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lstar = param.lstar;
a = param.a;
Gmult = param.Gmult;

lscaled = l/lstar;

Garray = Gmult*(1-a*lscaled).*exp(-lscaled);

end

