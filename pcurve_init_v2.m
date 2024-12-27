N = 100;
sfile = 'trillium2.mat';
sigma = linspace(0,1-1/N,N);

% Initial parameterization equally spaced in sigma 
z = (0:N-1)/N;
z = z';

x = zeros(N,1);
y = zeros(N,1);
dxs = zeros(N,1);
dys = zeros(N,1);

for j=1:N 
    [x(j), y(j), dxs(j), dys(j)] = cinit_xy(z(j));
end

  L = sum(sqrt((x(2:end)-x(1:end-1)).^2+ (y(2:end)-y(1:end-1)).^2));
  L = L + sqrt((x(1)-x(end)).^2+ (y(1)-y(end)).^2);

figure(1)
plot(x,y,'ko','LineWidth',2)
% for n = 1:numel(x)
%     text(x(n),y(n),num2str(n))
% end
axis equal 
title('Initial Parameterization')

tol = 1e-8;
it = 0;
residual = zeros(N-1,1);
% tres = sum(sqrt((x(2:end)-x(1:end-1)).^2+ (y(2:end)-y(1:end-1)).^2));
% tres = tres + sqrt((x(1)-x(end)).^2+ (y(1)-y(end)).^2);
% tres = tres - L;
% residual(1) = tres;
for j=1:N-1
    if j == N-1
        jp2 = 1;
    else
        jp2 = j+2;
    end
    residual(j) = (x(jp2) - x(j+1))^2 + (y(jp2)-y(j+1))^2 - ...
        (x(j+1)-x(j))^2 - (y(j+1)-y(j))^2;
end
res = max(abs(residual));

fprintf('Initial Newton residual %d \n',res);

if res < tol
    fprintf('Residual already less than tolerance \n')
end

% Save z and L in case of failure
z0 = z;
L0 = L;
    
while res > tol 
    % Try Newton 
    % Unknown order L then z2 to zN
    
    it = it+1;
      
    % Equations for z2 to zN
%     for j=1:N-1
%         if j == N-1
%             jp2 = 1;
%         else
%             jp2 = j+2;
%         end
%         residual(j) = (x(jp2) - x(j+1))^2 + (y(jp2)-y(j+1))^2 - ...
%             (x(j+1)-x(j))^2 - (y(j+1)-y(j))^2;
%     end
    A = zeros(N-1,N-1);
    for j=1:N-1
        if j == N-1
            jp2 = 1;
        else
            jp2 = j+2;
        end
        A(j,j) = 2*(x(j)-x(jp2))*dxs(j+1) + 2*(y(j)-y(jp2))*dys(j+1);
        if j ~= 1
            A(j,j-1) = +2*(x(j+1)-x(j))*dxs(j) + 2*(y(j+1)-y(j))*dys(j);
        end
        if j~= N-1
            A(j,j+1) = 2*(x(j+2)-x(j+1))*dxs(j+2) + 2*(y(j+2)-y(j+1))*dys(j+2);
        end
    end
    
    update = A\residual;
    L = L - update(1);
    z(2:end) = z(2:end) - update(1:end);
    
    for j=1:N 
        [x(j), y(j), dxs(j), dys(j)] = cinit_xy(z(j));
    end
    
    figure(2)
    plot(x,y,'ko','LineWidth',2)
    axis equal
    title('Newton Iteration')
%     for n = 1:numel(x)
%         text(x(n),y(n),num2str(n))
%     end
    getframe;
    % pause;
    
    for j=1:N-1
        if j == N-1
            jp2 = 1;
        else
            jp2 = j+2;
        end
        residual(j) = (x(jp2) - x(j+1))^2 + (y(jp2)-y(j+1))^2 - ...
            (x(j+1)-x(j))^2 - (y(j+1)-y(j))^2;
    end
    % residual
    res = max(abs(residual));
    fprintf('Newton residual %d \n',res);
    
    % Fail Newton step and do some local smoothing
    if res > 1 || it > 10
        fprintf('Failed Newton step, do some smoothing \n')
        it = 0;
        z = z0;
        L = L0;
        
        for j=1:N
            [x(j), y(j), dxs(j), dys(j)] = cinit_xy(z(j));
        end
        
        for j=1:N-1
            if j == N-1
                jp2 = 1;
            else
                jp2 = j+2;
            end
            residual(j) = (x(jp2) - x(j+1))^2 + (y(jp2)-y(j+1))^2 - ...
                (x(j+1)-x(j))^2 - (y(j+1)-y(j))^2;
        end
        res = max(abs(residual));
        fprintf('Initial smoothing Residual %d \n',res);
        
        % Smooth 50 times 
        for s=1:50
            % Smoothing means moving each point to the middle of its
            % neighbours
            for j=2:N
                stol = 1e-8;
                asres = 2*tol; 
                
                if j == N 
                    jn = 1;
                else 
                    jn = j+1;
                end
                zj = z(j);
                while asres>stol
                    sres = (x(jn)-x(j))^2 + (y(jn)-y(j))^2;
                    sres = sres - (x(j)-x(j-1))^2 - (y(j)-y(j-1))^2;
                    asres = abs(sres);
                    %fprintf('j= %d, residual %d\n',j,sres);
                    %pause
                    drive = -2*(x(jn)-x(j-1))*dxs(j);
                    drive = drive - 2*(y(jn)-y(j-1))*dys(j);
                    update = sres/drive;
                    zj = zj-update;
                    [x(j), y(j), dxs(j), dys(j)] = cinit_xy(zj);
                end
                % fprintf('j= %d iteration complete \n \n',j)
                z(j) = zj;
            end
            for j=1:N-1
                if j == N-1
                    jp2 = 1;
                else
                    jp2 = j+2;
                end
                residual(j) = (x(jp2) - x(j+1))^2 + (y(jp2)-y(j+1))^2 - ...
                    (x(j+1)-x(j))^2 - (y(j+1)-y(j))^2;
            end
            res = max(abs(residual));
            fprintf('Smoothing Residual %d \n',res);
        end
%         L = sum(sqrt((x(2:end)-x(1:end-1)).^2+ (y(2:end)-y(1:end-1)).^2));
%         L = L + sqrt((x(1)-x(end)).^2+ (y(1)-y(end)).^2);
               
        figure(3)
        plot(x,y,'ko','LineWidth',2)
        axis equal
        title('Smoothing Iteration')
%         for n = 1:numel(x)
%             text(x(n),y(n),num2str(n))
%         end
        getframe;
        % pause;
        % Save z and L in case of failure
        z0 = z;
        L0 = L;
    end    
end 

kappa = zeros(N,1);
kappass = zeros(N,1);

for j=1:N
    if j==1
        xp = x(N);
        yp = y(N);
    else 
        xp = x(j-1);
        yp = y(j-1);
    end
    
    if j==N
        xn = x(1);
        yn = y(1);
    else
        xn = x(j+1);
        yn = y(j+1);
    end
    xc = x(j);
    yc = y(j);
    
    kappa(j) = ((xn-2*xc+xp)*(yp-yn)+(yn-2*yc+yp)*(xn-xp))...
                /(2*L^3)*N^3;
end

for j=1:N
    if j==1
        kappap = kappa(N);
    else 
        kappap = kappa(j-1);
    end
    
    if j==N
        kappan = kappa(1);
    else
        kappan = kappa(j+1);
    end
    kappac = kappa(j);
    
    kappass(j) = (kappap-2*kappac+kappan)/L^2*N^2;
end

state0.x = x';
state0.y = y';
state0.L0 = L;
state0.L = L;
state0.kappa = kappa';
state0.kappass = kappass';

% save(sfile,'state0')

figure(4)
plot(sigma,kappa,'LineWidth',2)
title('Curvature')

figure(5)
plot(sigma,kappass,'LineWidth',2)
title('\kappa_{ss}')

fprintf('\n Curve initialized!\n')
