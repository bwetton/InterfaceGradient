function [GV, GVx, GVy, GVL, GVkappa] = Gcurve(L,x,y,kappa,param)

N = max(size(x));

% previous grid points
xp(2:N) = x(1:N-1);
xp(1) = x(N);
yp(2:N) = y(1:N-1);
yp(1) = y(N);

% next grid points
xn(1:N-1) = x(2:N);
xn(N) = x(1);
yn(1:N-1) = y(2:N);
yn(N) = y(1);

dsq = zeros(N,N);
dx = zeros(N,N);
dy = zeros(N,N);
ndx = zeros(N,N);
ndy = zeros(N,N);

for j=1:N
    for k = 1:N
        dx(j,k) = x(j)-x(k);
        dy(j,k) = y(j)-y(k);
        ndy(j,k) = yn(j) - yp(j);
        ndx(j,k) = xn(j) - xp(j);
        dsq(j,k) = (x(j)-x(k))^2 + (y(j)-y(k))^2;    
    end
end
alpha = -dx.*ndy+dy.*ndx;

G = Gfun(dsq,param);
Gprime = Gpfun(dsq,param);
Gpp = Gppfun(dsq,param);

%----------------------------------------------
%V1

% V1 row vector
Grow = sum(G');
V1 = 2*L/N*Grow;
V1 = V1.*kappa;

% Derivatives with L
GVL = V1/L;

% Derivatives with kappa
GVkappa = 2*L/N*Grow;

% Derivatives wrt x

J1x = diag(4*L*kappa/N.*sum(Gprime.*dx'));
for j=1:N
    for k=1:N
        if (j ~= k)
            J1x(j,k) = -4*L/N*kappa(j)*Gprime(j,k)*dx(j,k);
        end
    end
end

% Derivatives wrt y
J1y = diag(4*L*kappa/N.*sum(Gprime.*dy'));
for j=1:N
    for k=1:N
        if (j ~= k)
            J1y(j,k) = -4*L/N*kappa(j)*Gprime(j,k)*dy(j,k);
        end
    end
end

%----------------------------------------------
% V2

% V2 row vector
V2matrix = Gprime.*alpha;
V2 = -2*sum(V2matrix');

% Derivatives wrt x

J2diag = -4*Gpp.*dx.*(-dx.*ndy+dy.*ndx);
J2diag = J2diag + 2*Gprime.*ndy;
J2diag = sum(J2diag');
for j=1:N
    J2diag(j) = J2diag(j) - 2*Gprime(j,j)*ndy(j,j);
end

J2x = diag(J2diag);
Gdy = sum((Gprime.*dy)');

for j=1:N
    for k=1:N
        if (j~= k) 
            J2x(j,k) = 4*Gpp(j,k)*dx(j,k)*alpha(j,k) - 2*Gprime(j,k)*ndy(j,k);
        end
    end
    jn = j+1;
    if jn > N
        jn = 1;
    end
    J2x(j,jn) = J2x(j,jn) -2*Gdy(j);
    jp = j-1;
    if jp <1
        jp = N;
    end
    J2x(j,jp) = J2x(j,jp) + 2*Gdy(j);
end

% Derivatives with y

J2diag = -4*Gpp.*dy.*(-dx.*ndy+dy.*ndx);
J2diag = J2diag - 2*Gprime.*ndx;
J2diag = sum(J2diag');
for j=1:N
    J2diag(j) = J2diag(j) + 2*Gprime(j,j)*ndx(j,j);
end

J2y = diag(J2diag);
Gdx = sum((Gprime.*dx)');

for j=1:N
    for k=1:N
        if (j~= k) 
            J2y(j,k) = 4*Gpp(j,k)*dy(j,k)*alpha(j,k) + 2*Gprime(j,k)*ndx(j,k);
        end
    end
    jn = j+1;
    if jn > N
        jn = 1;
    end
    J2y(j,jn) = J2y(j,jn) +2*Gdx(j);
    jp = j-1;
    if jp <1
        jp = N;
    end
    J2y(j,jp) = J2y(j,jp) - 2*Gdx(j);
end

GV = V1+V2;
GVx = J1x+J2x;
GVy = J1y+J2y;

end

