function state_out = curve(state_in,finalT)
%
% Semi-group form of the curve evolution simulator
% Starts with input states, outpus states at finalT
% This version 2.0c specialized for the adhesion problem
%
% Computes the evolution of a curve undergoing geometric evolution.
% The normal velocity is user given:
%
% V = f(L,L0,kappa,kappa_{ss}) + [beta]
%
% The function f is found in the file curvef.m
% Additional adhesion terms are found in Gcurve, with the adhesion
%   function or distance and its derivatives in Gfun,Gpfun, Gppfun
% There is an option to add an area preserving term beta. This is set 
%    with the parameter "conserving" in curve_parameters.m
%
% Written by Brian Wetton, wetton@math.ubc.ca
% 
% v0.1 begun September 15, 2009.
% v1.0 begun May 31, 2011 - complete
%
% v1.1 begun June 5, 2011 - complete
% Added:
%    + proper error estimation
%    + option for periodic curves
%
% v2.0 begun May 28, 2024 
% Added: 
%    + Previouly Unrecorded: added adhesion term 
%       in Summer, 2023 (finished, but removed in this code) 
%    + Removed curve initialization from main code, 
%       allowing different approaches in separate codes.
%    + Changed code to semi-group form
%    + Removed periodic option
%    + Re-coded the aribtrary constant in the tangential 
%       velocity to be in integral form (in progress) 
%
% In the next versions, the following could be considered:
%    + after user trials, the code could be made more robust 
%    + add more options for curve initialization 
%

%% Initialization

% Load in parameters for the model and computations
param = curve_parameters;
sigma = param.sigma; % Local time step error tolerance

if param.movie == 1
    vid = VideoWriter(param.movie_name);
    open(vid);
end
state0 = state_in;

% Plot initial curve shape
figure(1)
xplot = state0.x;
yplot = state0.y;
plot(xplot, yplot, 'k-','LineWidth',3)
hold on
plot([xplot(1) xplot(end)],[yplot(1) yplot(end)], ...
    'k-', 'LineWidth',3)
hold off
axis equal
title('Initial Curve')

% Parameter array for later plotting
Ntot = max(size(state0.x));
z = (1:Ntot)/Ntot;

% Initialize the normal velocity if necessary 
if ~isfield(state0,'Vn')
    L = state0.L;
    L0 = state0.L;
    x = state0.x;
    y = state0.y;
    kappa = state0.kappa;
    kappass = state0.kappass;
    Ntot = max(size(x));
    Vn = zeros(1,Ntot);
    
    for i=1:Ntot
        [f, dL, dk, dkss] = curvef(L,L0,kappa(i),kappass(i),param);
        Vn(i) = f;
    end
    % Add G contribution to the normal velocity
    [GV GVx GVy GVL GVkappa] = Gcurve(L,x,y,kappa, param);
    Vn = Vn + GV; 
    state0.Vn = Vn;
end

if param.conserving == 1 && ~isfield(state0,'beta')
    Vn = state0.Vn;
    beta = sum(Vn)/Ntot;
    state0.beta = beta;
end

if isfield(state0,'k')
    k=state0.k;
else
    k = param.firstk;
end
done = 0;
fail=0;
success = 0;
t=0;

tv = t;
length = state0.L;
E = Energy(state0,param);
state2 = state0;

%% Time Stepping
% Perform time stepping to final time
fprintf('*** Start time stepping ***\n')
while done == 0
    % save old state in case of a failed time step
    state0 = state2;
    t0 = t;
    
    t1 = t + k;
    fail = 0;
    
    % First step
    % check for final time
    if t1> finalT 
        kold = k;
        t1 = finalT;
        t = t1;
        k = t1 - t0;
        if k < param.mink
            k = 2*param.mink;
        end
        done = 1;
    end
    % No prediction, just start with previous time step
    pstate = state0;
    
    % take a BE step
    state0.fail = 0;
    [state1, Nit1] = bestep(state0, k, param, pstate);
    if state1.fail == 1
        fail =1;
    end
    
    if done == 1
        state2 = state1;
    end
    
    % Second step
    if done ~=1 && fail == 0 
        t = t+k;
        % check for final time
        if t> finalT
            kold = k;
            t = finalT;
            k = t - t1;
            if k < param.mink
                k = 2*param.mink;
            end
            done = 1;
        end
        pstate = state1;
        % take a BE step
        state1.fail = 0;
        [state2, Nit2] = bestep(state1, k, param, pstate);
        if state2.fail == 1
            fail =1;
        end
    end
    
    % Coarse step for error estimation
    if done ~=1 && fail == 0
        % do a prediction step
        %pstate = fepredictor(state0,2*k,param);
        pstate = state2;
        % take a BE step
        state0.fail = 0;
        [state3, Nit3] = bestep(state0, 2*k, param, pstate);
        if state3.fail == 1
            fail =1;
        end
        % Estimate of the error
        dx = state2.x-state3.x;
        dy = state2.y-state3.y;
        dkappa = state2.kappa-state3.kappa;
        dkappass = state2.kappass-state3.kappass;
        Eest = max(abs([dx dy dkappa dkappass]));
        Eest = max(abs([dx dy]));
        
        if Eest>sigma
            fprintf('Failing time step loss of accuracy \n')
            fail = 1;
        end
    end
    
    %% Time step done, check fail or success 
    if fail == 1
        t = t0;
        state = state0;
        k = k/2;
        success = 0;
        done =0;
        fprintf('BE step fail, reducing step size to %d\n', k)
    else 
        fprintf('Time step %10.3e, (%d iterations) \n',k,Nit1+Nit2+Nit3)
        tv = [tv t];
        length = [length state2.L];
        E = [E Energy(state2,param)];

        kadjust = min(param.alpha,0.8*sqrt(sigma/Eest));
        k = k*kadjust;
        if k > param.maxk
            k = param.maxk;
        end
        
        % Plot curve shape
        figure(2)
        xplot = state2.x;
        yplot = state2.y;
        plot(xplot, yplot, 'k-','LineWidth',3)
        hold on
        plot([xplot(1) xplot(end)],[yplot(1) yplot(end)], ...
            'k-', 'LineWidth',3)
        hold off
        axis equal
        title(sprintf('Time %d',t))
        
%         %Keith adds tracking center of mass of region
%         figure(3)
%         sx=size(xplot);
%         xnum=sx(2);
%         st=size(tv);
%         knum=st(2);
%         if knum<100
%           clear avx
%           clear avy
%         end
%         avx(knum)=sum(xplot)/xnum;
%         avy(knum)=sum(yplot)/xnum;
%         rad=sum(sqrt((xplot-avx(knum)).^2+(yplot-avy(knum)).^2))/xnum;
%         pt=(0:0.05:2*pi);
%         % End Keith changes 
%         plot(xplot, yplot, 'k-','LineWidth',3)
%         hold on
%         plot([xplot(1) xplot(end)],[yplot(1) yplot(end)], ...
%             'k-', 'LineWidth',3)
%         plot(avx,avy,'b','linewidth',2)
%         plot(avx(end),avy(end),'ro','markersize',15)
%         plot(rad*cos(pt)+avx(st(2)),rad*sin(pt)+avy(st(2)),'b.')
%         hold off
%         axis equal
%         title(sprintf('Time %d',t))
        frame = getframe(gcf);
        if param.movie == 1
            writeVideo(vid,frame);
        end
    end
    if k < param.mink 
        fprintf('*** Time step below minimum - solution fail ***\n')
        done = 1;
        fail =1;
    end
    
    figure(4)
    %rho=param.rho;
    L0=state2.L0;
    plot(tv, length/L0, 'LineWidth',2)
    title('L/L_0 vs time')
    
    figure(6)
    plot(tv, E, 'LineWidth',2)
    title('Energy vs time')
    
    figure(5)
    kappa=state2.kappa;
    plot(z,kappa,'bo','markersize',5)
    hold on
    plot(z,kappa,'b')
    hold off
    title(sprintf('Curvature at Time %d',t))
end

%% Wrap up final state

state2.k = kold;
state_out = state2;

if fail == 0 
    fprintf('*** Time stepping complete ***\n')
    fprintf('*** %d total time steps \n',max(size(tv))-1)
end
if param.movie == 1
    close(vid);
end


