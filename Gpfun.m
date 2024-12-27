function Gparray = Gpfun(l,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lstar = param.lstar;
a = param.a;
Gmult = param.Gmult;

lscaled = l/lstar;

Gparray = Gmult*(-a/lstar+a*lscaled/lstar-1/lstar).*exp(-lscaled);

end