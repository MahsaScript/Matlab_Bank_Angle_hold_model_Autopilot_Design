%% Bank Angle hold model of  Autopilot Design

% Clear data and figures
clc
clear
close all

% Define State Space Data
A=[0 500 0 -500 0 0 0;
   0 0 0 0 1 0 0;
   0.000108 -32.17 -0.01328 -7.326 -1.196 0.001565 0.07397;
   2.076e-06 -7.394e-13 0.0002553 -0.6398 0.9378 -2.445e-07 -0.001357;
   2.83e-20 0 -3.48e-18 -1.568 -0.8791 0 -0.1137; 
   0 0 0 0 0 -1 0; 
   0 0 0 0 0 0 -20.2];
B=[0 0
   0 0
   0 0
   0 0
   0 0
   1 0
   0 20.2 ];
C= [1 0 0 0 0 0 0;
    0 57.3 0 0 0 0 0;
    0 0 1 0 0 0 0;
    0 0 0 57.3 0 0 0;
    0 0 0 0 57.3 0 0];
D=[0 0;
   0 0;
   0 0;
   0 0;
   0 0];

trim_control_lo=[-2.4607;
                 0;
                 0];
trim_thrust_lo=2.1206e+03;
Trim_state_lo=1e4*[0 0 1.5 0 0 0 0.05 0 0 0];

%% (a) Open-Loop System

G=ss(A,B,C,D,'statename',{'h','theta','v','alpha','q','delta_t','delta_e'});
display(G)

%% (b) Open-Loop Analysis

%Stability analysis
lambda = eig(A);
damp(lambda)
% Pole and Zero Plot
figure
pzplot(G)

% Step Response of System
figure
step(G,0.5)
% Impulse Response of System
figure
impulse(G,0.5)
% Bode Plot
figure
bode(G)
% Controllability and Observability
Co = rank(ctrb(G));
if Co==7
    fprintf('\nSystem is Controllable\n')
end
Ob = rank(obsv(G));
if Ob==7
    fprintf('System is Observable\n\n')
end

% Simulink
tmax = 0.5; % Simulation stop time
SYS = 'my_plant'; % Model’s name
so=sim(SYS); % Run the simulation

% Plot Simulation Results
Dataout=so.simout;
Datain=so.simin;
figure
plot(so.tout,Datain,so.tout,Dataout,'linewidth',2)
grid on
xlabel('s')
ylabel('u and y')
title('Step Response Simulink')

%% (c) State-Feedback Tracking Controller

% State-feedback tracking controllers Using LQR
Q = C'*C;
R = 1;
% lqr controller
K = lqr(A,B,Q,R);
Nbar = rscale(A,B,C(1:2,:),[0 0;0 0],K)/1.43;
% Initialize State Matrix
Anw = (A-B*K);
Bnw = B;
Cnw = C;
Dnw = D;
% state space model
G_cl=ss(Anw,Bnw.*Nbar',Cnw,Dnw);
% time vector
t = 0:0.1:500;
% Input response
r =[15000*ones(numel(t),1) zeros(numel(t),1)];
% Simulate model
[y,t,x]=lsim(G_cl,r,t);
% Plot result
figure
plot(t,y(:,1))
grid on
ylabel('Height (ft)')
xlabel('Time (seconds)')
title('Height Vs Time')
info=stepinfo(y(:,1),t);
title(['Using LQR, OverShoot:',num2str(info.Overshoot),' %'])

% State-feedback tracking controllers Using Place
% selected Pole 
p1 = -18;
p2 = -9;
p3 = -1.0129 + 1.3388i;
p4 = -1.0129 - 1.3388i;
p5 = -1.5;
p6 = -1+0.5i;
p7 = -1-0.5i;
% place pole
K = place(A,B,[p1 p2 p3 p4 p5 p6 p7]/2);
Nbar = rscale(A,B,C(1:2,:),[0 0;0 0],K)/(1.3*10^4);
% Initialize State Matrix
Anw = (A-B*K);
Bnw = B;
Cnw = C;
Dnw = D;
% state space model
G_cl=ss(Anw,Bnw.*Nbar',Cnw,Dnw);
% time vector
t = 0:0.1:100;
% Input response
r =[15000*ones(numel(t),1) zeros(numel(t),1)];
% Simulate Model
[y,t,x]=lsim(G_cl,r,t);
% Plot result
figure
plot(t,y(:,1))
grid on
ylabel('Height (ft)')
xlabel('Time (seconds)')
title('Height Vs Time')
info=stepinfo(y(:,1),t);
title(['Using Place, OverShoot:',num2str(info.Overshoot),' %'])

%% (d) Close-Loop Performance

% statespace to transfer function model
[b,a]=ss2tf(G.A,G.B,G.C,G.D,2);
% define tf model
sys1=tf(b(1,:),a);
sys2=tf(b(2,:),a);
sys3=tf(b(3,:),a);
sys4=tf(b(4,:),a);
sys5=tf(b(5,:),a);
% Transfer function Model
sys=[sys1;sys2;sys3;sys4;sys5];
% Feedback controller
s=tf('s');
K=-0.1e-04*((1+10*s)*(1+200*s))/s;
% Overall system
sys_cl=feedback(K*sys,1,1,1,-1);
% time vector
t = 0:0.5:500;
% Input response
r =15000*ones(numel(t),1);
% Simulate Model
[y,t,x]=lsim(sys_cl,r,t);
% Plot result
figure
plot(t,y(:,1))
grid on
ylabel('Height (ft)')
xlabel('Time (seconds)')
title('Height Vs Time')

%Stability analysis
damp(sys_cl)
% Pole and Zero Plot
figure
pzplot(sys_cl)
% Bode Plot
figure
bode(sys_cl)


[n1, nu] = size(B);n2 = 1; % One integrator added
E = [0 1 0 0 0 0 0 0]; % E such that E*x = the state you want to track
% Augmented system
Aaug = [A zeros(n1, n2);-E]; 
Baug = [B;zeros(n2,nu)];
% Design state-feedback for the augmented system (Aaug, Baug)
K = lqr(A,B,Q,R);
% controller gain
K=[K,[0;0],[0;0]];
% Simulate Model
data=sim('Plant_Close.slx');
% Plot result
figure
plot(data.tout,data.y.Data(:,1))
grid on
ylabel('Height (ft)')
xlabel('Time (seconds)')
title('Using Integral Controller')


%  variation from nomial value of system
Aaug=1.2*Aaug;
Baug=1.2*Baug;
K=1.2*K;
% Simulate Model
data=sim('Plant_Close.slx');
% Plot result
figure
subplot(121)
plot(data.tout,data.y.Data(:,1))
grid on
ylabel('Height (ft)')
xlabel('Time (seconds)')
title('Height vs time')
subplot(122)
plot(data.tout,data.y.Data(:,3))
grid on
ylabel('Velocity (ft/s)')
xlabel('Time (seconds)')
title('Velocity vs time')
sgtitle('For 20% Variation')

%% (g) Controller in Discrete Time

% Sampling Period
Td=0.2;
% Integral Controller Set-up
[n1, nu] = size(B);n2 = 1; % One integrator added
E = [0 1 0 0 0 0 0 0]; % E such that E*x = the state you want to track
% Augmented system
Aaug = [A zeros(n1, n2);-E]; 
Baug = [B;zeros(n2,nu)];
% Design state-feedback for the augmented system (Aaug, Baug)
K = lqr(A,B,Q,R);
% controller gain
K=[K,[0;0],[0;0]];
% Simulate Model
data=sim('Plant_discrete');
% Plot result
figure
plot(data.tout,data.y.Data(:,1))
grid on
ylabel('Height (ft)')
xlabel('Time (seconds)')

function[Nbar]=rscale(a,b,c,d,k)
A=a; B=b; C=c; D=d; K=k;
% compute Nbar
s = size(A,1);
Z = [zeros([1,s]) 1 1];
N = inv([A,B;C,D])*Z';
Nx = N(1:s);
Nu = N(1+s);
Nbar=Nu + K*Nx;
end