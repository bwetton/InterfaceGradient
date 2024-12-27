function [state, iter] = bestep(state_in, k, param, pstate)

state = state_in;

L = pstate.L;
x = pstate.x;
y = pstate.y;
xold = state.x;
yold = state.y;
kappa = pstate.kappa;
kappass = pstate.kappass;
Vn = pstate.Vn;
if param.conserving == 1
    beta = state.beta;
end

Ntot = max(size(x));


% The size of the linear system for Newton's method
Nsolve = 5*Ntot +1;
% The order of unknowns in this system are L then x, y, kappa, kappass,
% and Vn at each successive point.

tolerance = param.tol;
residual = 2*tolerance;
fail = 0;
iter = 0;


while (residual > tolerance) && (fail == 0) 
    [GV, GVx, GVy, GVL, GVkappa] = Gcurve(L,x,y,kappa,param);
    
    if param.conserving == 0
        A = zeros(Nsolve, Nsolve);
        res = zeros(Nsolve,1);
    else
        A = zeros(Nsolve+1, Nsolve+1);
        res = zeros(Nsolve+1,1);
    end
    dx = L/Ntot;
 %% Preliminaries to handle the periodicity    
    for i=1:Ntot
        if i==1 
            xp = x(Ntot);
            yp = y(Ntot);
            kappap = kappa(Ntot);
            ip = 5*Ntot-3;
        else
            xp = x(i-1);
            yp = y(i-1);
            kappap = kappa(i-1);
            ip = 5*(i-1)-3;
        end
        
        if i==Ntot
            xn = x(1);
            yn = y(1);
            kappan = kappa(1);
            in = 2;
        else
            in = 5*(i+1)-3;
            xn = x(i+1);
            yn = y(i+1);
            kappan = kappa(i+1);
        end
        
        ic = 5*i-3;
        xc = x(i);
        yc = y(i);
        kappac = kappa(i);

%% length residual for y equations 
        length2 = (xc-xp)^2 + (yc-yp)^2;
        res(ic+1) = length2-L^2/Ntot^2;
        A(ic+1,ip) = -2*(xc-xp);
        A(ic+1,ip+1) = -2*(yc-yp);
        A(ic+1,ic) = 2*(xc-xp);
        A(ic+1,ic+1) = 2*(yc-yp);
        A(ic+1,1) = -2*L/Ntot^2;
               
        dy1 = yc-yp;
        dy2 = yn-yc;
        dyc = dy2+dy1;
        
%% BE for normal motion for x equations 
        res(ic) = ((xc - xold(i))*(-dyc) + (yc - yold(i))*(xn-xp)) ...
            - 2*dx*Vn(i)*k;
        A(ic, 1) = -2*Vn(i)*k/Ntot;
        A(ic, ip) = -(yc - yold(i));
        A(ic, ip+1) = (xc - xold(i));
        A(ic, ic) = (-dyc);
        A(ic, ic+1) = (xn-xp);
        A(ic, in) = (yc - yold(i));
        A(ic, in+1) = -(xc - xold(i));
        A(ic, ic+4) = -2*dx*k;
        
        if param.conserving == 1
            res(ic) = res(ic)+2*dx*beta*k;
            A(ic, 1) = A(ic,1) + 2*beta*k/Ntot;
            A(ic, Nsolve+1) = 2*dx*k;
        end
       
%% Fix sum of W=0 for L equation
        
        res(1) = res(1) + (xc - xold(i))*(xn-xp) + (yc-yold(i))*(dyc);
        A(1, ip) = A(1,ip) -(xc - xold(i));
        A(1, ip+1) = A(1,ip+1)-(yc - yold(i));
        A(1, ic) = A(1,ic) +(xn-xp);
        A(1, ic+1) = A(1,ic+1) + (dyc);
        A(1, in) = A(1,in) + (xc - xold(i));
        A(1, in+1) = A(1,in+1)+(yc-yold(i));
        
%% kappa
        res(ic+2) = ((xp+xn-2*xc)*(-dyc)/2 + (dy2-dy1)*(xn-xp)/2)/dx^3 - ...
            kappa(i);
        A(ic+2, ic+2) = -1;
        A(ic+2, 1) = -3*((xp+xn-2*xc)*(-dyc)/2 + (dy2-dy1)*(xn-xp)/2)* ...
            Ntot^3/L^4;
        A(ic+2, in) = (-2*dy1)/2/dx^3;
        A(ic+2, in+1) = (-(xp+xn-2*xc)+(xn-xp))/2/dx^3;
        A(ic+2, ic) = (dyc)/dx^3;
        A(ic+2, ic+1) = (-(xn-xp))/dx^3;
        A(ic+2, ip) = ((-2*dy2))/2/dx^3;
        A(ic+2, ip+1) = ((xp+xn-2*xc)+(xn-xp))/2/dx^3;
        
%% kappass
        res(ic+3) = (kappap+kappan-2*kappac)/dx^2- kappass(i);
        
        A(ic+3,ic+3) = -1;
        A(ic+3,ip+2) = 1/dx^2;
        A(ic+3,in+2) = 1/dx^2;
        A(ic+3, ic+2) = -2/dx^2;
        A(ic+3,1) = -2*(kappap+kappan-2*kappac)*Ntot^2/L^3;
        
%% normal velocity
        
        [f, df1, df2, df3] = curvef(L,state.L0,kappa(i), kappass(i),param);
        f = f+GV(i);
        res(ic+4) = f - Vn(i);
        
        A(ic+4, 1) = df1+GVL(i);
        A(ic+4, ic+2) = df2+GVkappa(i);
        A(ic+4, ic+3) = df3;
        A(ic+4, ic+4) = -1;
        
        % add in velocity jacobian entries for x,y from G term 
        for j=1:Ntot
            A(ic+4,5*j-3) = GVx(i,j);
            A(ic+4,5*j-2) = GVy(i,j);
        end
        
    end
    
    if param.conserving == 1
        res(Nsolve+1) = sum(Vn)/Ntot-beta;
        for i=1:Ntot
            A(Nsolve+1, 5*i+1) = 1/Ntot;
        end
        A(Nsolve+1,Nsolve+1) = -1;
    end
    
    residual = max(abs(res));
    % fprintf('Residual: %d\n', residual)
    
    update = (A\res)';
    
    L = L - update(1);
    x = x - update(2:5:(Nsolve-4));
    y = y - update(3:5:(Nsolve-3));
    kappa = kappa - update(4:5:(Nsolve-2));
    kappass = kappass - update(5:5:(Nsolve-1));
    Vn = Vn - update(6:5:Nsolve);
    if param.conserving == 1
        beta = beta - update(Nsolve+1);
    end
    
    iter = iter +1;
    if iter > param.maxNewton 
        fprintf('Failing time step, maximum Newton iterations exceeded \n')
        fail =1;
    end
    
    msize = max(max(abs(x)),max(abs(y)));
    if msize > param.maxSize 
        fprintf('Failing time step, wild Newton iterate \n')
        fail = 1;
    end
end

if fail == 1
    state.fail = 1;
else
    state.fail = 0;
    state.x = x;
    state.y = y;
    state.kappa = kappa;
    state.kappass = kappass;
    state.Vn = Vn;
    state.L = L;
    if param.conserving == 1
        state.beta = beta;
    end
end

% fprintf('Newton iteration successful with %d iterations \n',iter)

return
