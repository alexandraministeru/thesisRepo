function [u,y,sys_d] = getIOData(Ts)
%getIOData(Ts) Function used to define a mass-spring-damper 2nd order system
%and get IO open-loop data sampled with sampling period Ts.
%
% Input arguments
%----------------
% Ts    : Sampling period (in seconds).
%
% Output arguments: 
%------------------
% u     : input vector - external force acting on the body (in Newtons)
% y     : output vector - displacement of the body (in meters)
% sys_d : discretized system
%==========================================================================

% System parameters
m = 1; % mass
k = 10; % spring constant
b = 0.8; % damping constant

% State space representation
A_sys = [0 1; -k/m -b/m];
B_sys = [0; 1/m];
C_sys = [1 0];
D_sys = 0;
sys = ss(A_sys,B_sys,C_sys,D_sys);

% Discretized system
sys_d = c2d(sys,Ts);

% Time vector
t = 0:Ts:20;

% PRBS input for persistency of excitation
u = idinput(length(t),'PRBS',[0 1],[0 2]);

% Simulate system and get open-loop I/O data set
y = lsim(sys_d,u,t); % zero initial condition
figure
plot(t,y);
grid on
xlabel('Time (in s)')
ylabel('Displacement')
title('Open-loop response to PRBS input')
set(gcf,'Color','White')

end