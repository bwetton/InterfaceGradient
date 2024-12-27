function E = Energy(state2,param)

kappa=state2.kappa;
L0=state2.L0;
eps=param.epsilon;
mu=param.mu;
rho=param.rho;
L=state2.L; 
x=state2.x;
y=state2.y;
Ntot = max(size(kappa));
Gterm=0;
for kk=1:Ntot
 xp=x-x(kk);
 yp=y-y(kk);
 %ell=sqrt(xp.^2+yp.^2);
 ell=xp.^2+yp.^2;
 Gs=Gfun(ell,param);
 Gterm=Gterm+sum(Gs)*L/Ntot;
end
Gterm=Gterm*L/Ntot;
%[W,Wp,Wpp,Wppp,Wpppp]=Wfun(kappa,param);
E=Gterm+sum(eps*kappa.^2/2)*L/Ntot+mu/2*(L-rho*L0)^2;
