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

%% Set indexes
j = length(y);
M = 10;
N = 10;
Nbar = j-M-N+1;
i = 1;
m = 1;
l = 1;
%% Data matrices
% Past data
Up = constructHankelMat(u,i,M,Nbar);
Yp = constructHankelMat(y,i,M,Nbar);
Wp = [Yp;Up];

% Future data
Uf = constructHankelMat(u,i+M,N,Nbar);
Yf = constructHankelMat(y,i+M,N,Nbar);

%% LQ decomposition
H = [Wp; Uf; Yf];
R = triu(qr(H',0));
R = R';

R11 = R(1:2*M, 1:2*M);
R21 = R(2*M+1:2*M+N, 1:2*M);
R22 = R(2*M+1:2*M+N, 2*M+1:2*M+N);
R31 = R(2*M+N+1:end, 1:2*M);
R32 = R(2*M+N+1:end, 2*M+1:2*M+N);
R33 = R(2*M+N+1:end, 2*M+N+1:end);

%% Solve least squares
L = [R31 R32]/[R11 zeros(2*M,M); R21 R22];
Lw = L(:,1:N*(m+l));
Lu = L(:,N*(m+l)+1:end);

%% Take SVD and approximate Lw
[U_svd,S_svd,V_svd] = svd(Lw);
figure(1)
sv = diag(S_svd);
% semilogy(sv,'x')
plot(sv,'x')
title('Singular values')
grid on
[n,~] = ginput(1);close
n = round(n); 
Lw = U_svd(:,1:n)*S_svd(1:n,1:n)*V_svd(1:n,:);
 
%% Control loop
tFinal = 7;
x0 = [1 0.2]; % initial condition

% Define weights
Q_spc = eye(N);
R_spc = eye(N);

% Keep track of states
x = zeros(2,tFinal+1);
x(:,1) = x0;

% Keep track of outputs
out = zeros(tFinal,1);

% Latest M known inputs
yp = y(j-M+1:end);
up = u(j-M+1:end);
wp = [yp; up];

for t=1:tFinal
    % Control law
    uf = (R_spc + Lu'*Q_spc*Lu)\(Lu'*Q_spc*(Lw*wp));
    uStar = uf(1);

    % Simulate output
    x(:,t+1) = A_sys*x(:,t) + B_sys.*uStar;
    out(t) = C_sys*x(:,t) + D_sys.*uStar;

    % Update IO data with the most recent values
    yp = [yp(2:end); out(t)];
    up = [up(2:end); uStar];
    wp = [yp; up];
end

figure
plot(1:tFinal,out)
xlabel('Time(s)')
ylabel('Displacement')
title('Subspace predictive control')
grid on


