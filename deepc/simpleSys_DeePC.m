clearvars;clc;
%% Define system, get I/O data
% Mass-spring-damper 2nd order system
m = 1;
k = 5;
b = 0.02;

% State space representation
A_sys = [0 1; -k/m -b/m];
B_sys = [0; 1/m];
C_sys = [1 0];
D_sys = 0;
sys = ss(A_sys,B_sys,C_sys,D_sys);

% Time vector
t = 0:0.1:200;
% PRBS input for persistency of excitation
u = idinput(length(t),'PRBS',[0 1],[0 1]); 

% Simulate system and get open-loop I/O data set
y = lsim(sys,u,t);
figure
plot(t,y);
grid on
xlabel('Time(s)')
ylabel('Displacement')
title('Open-loop response to PRBS input')

%% DeePC parameters
N = 100;
p = 10;
f = 10;
Nbar = N-p-f+1;
i = 1;

%% Data matrices
% Past data
U_ipNbar = constructHankelMat(u,i,p,Nbar);
Y_ipNbar = constructHankelMat(y,i,p,Nbar);

% Future data
U_ipfNbar = constructHankelMat(u,i+p,f,Nbar);
Y_ipfNbar = constructHankelMat(y,i+p,f,Nbar);

%% Solve the constrained optimization problem
% yalmiprootdirectory = 'D:\lib\YALMIP';
% addpath(genpath(yalmiprootdirectory));
% addpath(genpath('D:\Program Files\MOSEK\10.0\toolbox\r2017aom'));

% Define weights
Q = eye(f);
R = eye(f);

% Past data for prediction
U_ihatp1 = constructHankelMat(u,i+N-p,p,1);
Y_ihatp1 = constructHankelMat(y,i+N-p,p,1);

tFinal = 7;
x0 = [1 0.2]; % initial condition

% Keep track of states
x = zeros(2,tFinal+1);
x(:,1) = x0;

% Keep track of output
out = zeros(tFinal,1);

% Control loop
for t=1:tFinal
    % % Define decision variables
    % U_iphatf1 = sdpvar(f,1);
    % g = sdpvar(Nbar,1);
    % 
    % % Define constraints
    % A = [U_ipNbar;
    %     Y_ipNbar;
    %     U_ipfNbar];
    % 
    % b = [U_ihatp1;
    %     Y_ihatp1;
    %     U_iphatf1];
    % 
    % con = [A*g==b];
    % 
    % % Define objective function - zero reference tracking
    % obj = (Y_ipfNbar*g)'*Q*(Y_ipfNbar*g) + U_iphatf1'*R*U_iphatf1;
    % 
    % % Options
    % options = sdpsettings('verbose',1); % implicitly uses quadprog, non-convex problem(?), thus no optimization
    % 
    % % Solve optimization problem
    % sol = optimize(con,obj,options);
    % 
    % % Assign the solution    
    % if sol.problem == 0
    %     U_iphatf1 = value(U_iphatf1);
    %     g = value(g);
    % else
    %     disp('Something went wrong');
    %     sol.info
    %     yalmiperror(sol.problem)
    % end
    % 
    % % Apply first input
    % uStar = U_iphatf1(1);
    
    % Decision variables
    g = zeros(Nbar,1);
    U_iphatf1 = zeros(f,1);
    Y_iphatf1 = zeros(f,1);   
    z = [g; U_iphatf1; Y_iphatf1];

    % Equality constraints
    Aeq = [U_ipNbar zeros(p) zeros(p);
        Y_ipNbar zeros(p) zeros(p);
        U_ipfNbar -eye(f) zeros(f);
        Y_ipfNbar zeros(f) -eye(f)];

    beq = [U_ihatp1;
        Y_ihatp1;
        zeros(f,1);
        zeros(f,1)];

    % Objective function
    objFcn = @(z) (z'*[zeros(Nbar) zeros(Nbar,f) zeros(Nbar,f);
                        zeros(f,Nbar) R zeros(f);
                        zeros(f,Nbar) zeros(f) Q]*z);

    options = optimoptions('fmincon','Algorithm','sqp');
    
    % Solve optimization problem
    [z,fval,exitflag,output] = fmincon(objFcn,z,[],[],Aeq,beq,[],[],[],options);

    g = z(1:Nbar);
    U_iphatf1 = z(Nbar+1:Nbar+f);
    Y_iphatf1 = z(Nbar+f+1:end);
    
    % Choose only first input
    uStar = U_iphatf1(1);
    
    % Apply optimal input, simulate output
    x(:,t+1) = A_sys*x(:,t) + B_sys.*uStar;
    out(t) = C_sys*x(:,t) + D_sys.*uStar;
    
    % Update past data with most recent I/O data
    U_ihatp1 = [U_ihatp1(2:end); uStar];
    Y_ihatp1 = [Y_ihatp1(2:end); out(t)];
end

figure
plot(1:tFinal,out)
xlabel('Time(s)')
ylabel('Displacement')
title('DeePC')
grid on

