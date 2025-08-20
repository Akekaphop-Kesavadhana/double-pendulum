close all; clear; clc; 

% DPIC - Double Pendulum Inverted on a Cart
% Akekaphop Kesavadhana
% 03/26/2024

% The objective of this MATLAB code was to
% reproduce the results of 'Control Theory: The Double Pendulum
% Inverted on a Cart' by Ian J. P. Crowe-Wright. Furthermore, 
% the motivation was to explore the plethora of control
% techniques applied to nonlinear dynamical systems and act as the beginning to 
% constructing an in-person double pendulum inverted on a cart bench apparatus for 
% further testing and application of linear control theory on chaotic nonlinear systems.

% Future objectives moving further in the project is to develop a swing up controller 
% that will utilize a Lyapunov function that defines the energy within the
% system and thus be able to allow a controller to swing both pendulums up
% afterwhich the controller mode will switch from swing up to stabilization
% and maintain the zero angle reference.

% Citation/Reference :
% Crowe-Wright, Ian J P. "Control Theory: The Double Pendulum Inverted on a Cart." (2018). 
% https://digitalrepository.unm.edu/math_etds/132

% Constant initializations and definitions
g  = 9.8;               % Gravitational constant in [m/s^2]
m0 = 1.5;               % Mass of cart in [kg]
m1 = 0.5;               % Mass of lower pendulum link in [kg]
m2 = 0.75;              % Mass of upper pendulum link in [kg]
L1 = 0.5;               % Length of lower pendulum link in [m]
L2 = 0.75;              % Length of upper pendulum link in [m]
l1 = (1/2)*L1;          % Distance between the base of first pend. link and its center of mass in [m]
l2 = (1/2)*L2;          % Distance between the base of second pend. link and its center of mass in [m]
I1 = (1/12)*m1*(L1)^2;  % Moment of inertia of first pendulum link in [kg*(m)^2]
I2 = (1/12)*m2*(L2)^2;  % Moment of inertia of second pendulum link in [kg*(m)^2]

% MATLAB Matrix standard [row x column]

% D*theta_double_dot + C*theta_dot + G = H*u -> Second order nonlinear system matrix model
% x_dot = A*x + B*u + L -> First order nonlinear system matrix model
% x_dot = A*x + B*u -> First order linearized system matrix model (Standard state space form)

% D = [m0+m1+m2                (m1*l1+m2*L2)*cos(theta1)     m2*l2*cos(theta2);
%     (m1*l1+m2*L1)*cos(theta1) m1*(l1)^2+m2*(L1)^2+I1       m2*L1*l2*cos(theta1-theta2);
%      m2*l2*cos(theta2)        m2*L1*l2*cos(theta1-theta2)  m2*(l2)^2+I2];

% C = [0     -(m1*l1+m2*L1)*sin(theta1)*theta1_dot    -m2*l2*sin(theta2)*theta2_dot;
%      0                      0                        m2*L1*l2*sin(theta1-theta2)*theta2_dot;
%      0     -m2*L1*l2*sin(theta1-theta2)*theta1_dot   0];

% G = [0;
%    -(m1*l1+m2*L1)*g*sin(theta1);
%     -m2*g*l2*sin(theta2)];

% H = [1;
%      0;
%      0];

% NOTE: D is symmetric and nonsingular -> D^-1 exists and is also symmetric

% D*theta_double_dot + C*theta_dot + G = H*u -> Convert second order to first order system (solve for theta_double_dot)
% D*theta_double_dot = -C*theta_dot - G + H*u
% theta_double_dot = -D^-1*C*theta_dot - D^-1*G + D^-1*H*u

% Let: x = [theta          Let: x_dot = [theta_dot
%           theta_dot]                   theta_double_dot]
% x_dot = A*x + B*u + L -> A,B, and L are functions of x

% A = [0    I;
%      0   -D^-1*C];
% 
% B = [0;
%      D^-1*H];
% 
% L = [0;
%     -D^-1*G];

% The above is a nonlinear problem of the form 
% x_dot(t) = f(x(t),u(t)), x(0) = x0

%  In order to reduce the system to a linear form, series expansions about
%  theta = 0 are employed, its is assumed that the pendulum is close to the
%  vertical. The L term is absorbed into the x coeffient. This eliminates
%  theta dependence. Technically, this is a linearization about x=0. The
%  linear system could also be found by computing the Jacobian and
%  evaluating at a fixed point.

% x_dot = A*x + B*u
% C = eye(6); % y = Cx, C = I (Identity matrix), thus, y = x, thus y = x which means the system is fully observable, ie each component of the state vector can be measured
% D = [0 0 0 0 0 0]';

%% ========= Linearized state space form and matrix constituents =========

p = 1/(4*m0*m1 + 3*m0*m2 + (m1)^2 + m1*m2);

a42 = -(3/2)*p*(2*(m1)^2 + 5*m1*m2 + 2*(m2)^2)*g;
a43 =  (3/2)*p*m1*m2*g;
a52 =  (3/2)*(p/L1)*(4*m0*m1 + 8*m0*m2 + 4*(m1)^2 + 9*m1*m2 + 2*(m2)^2)*g;
a53 = -(9/2)*(p/L1)*(2*m0*m2 + m1*m2)*g;
% a62 = -(9/2)*((p*g)/L2)*(2*m0*m1 + 4*m0*m2 + (m1)^2 + 2*m1*m2);
a62 = -(9*p*g*(2*m0+m1)*(m1+2*m2))/(2*L2); % missing g - gravitational constant, a multiplicative factor of 2*m2 is a result of common denominators
% a63 =  (3/2)*(p/L2)*((m1)^2 + 4*m0*m1 + 12*m0*m2 + 4*m1*m2); 
a63 =  (3/2)*((p*g)/L2)*((m1)^2 + 4*m0*m1 + 12*m0*m2 + 4*m1*m2); % missing g - gravitational constant

b4 = p*(4*m1 + 3*m2);
b5 =  -((3*p)/L1)*(2*m1 + m2);
b6 =   (3*p*m1)/L2;

   % x    x_dot  theta1 theta1_dot theta2 theta2_dot
A = [0    0      0      1          0      0;        % System matrix       s = [x
     0    0      0      0          1      0;                              %    x_dot
     0    0      0      0          0      1;                              %    theta1
     0    a42    a43    0          0      0;                              %    theta1_dot
     0    a52    a53    0          0      0;                              %    theta2
     0    a62    a63    0          0      0];                             %    theta2_dot]

B = [0;
     0;
     0;
     b4;
     b5;
     b6];

C = eye(6); 

D = [0 0 0 0 0 0]';


% THESIS RESULTS : 

% A = [0    0      0       1 0 0;
%      0    0      0       0 1 0;
%      0    0      0       0 0 1;
%      0   -7.35   0.7875  0 0 0;
%      0    73.5  -33.0750 0 0 0;
%      0   -58.8   51.1    0 0 0];
% 
% B = [0;
%      0;
%      0;
%      0.6071;
%     -1.5;
%      0.2857];

% MY ORIGINAL RESULTS (*EQNS GIVEN WERE INCORRECT) : 

% A =      0         0         0    1.0000         0         0
%          0         0         0         0    1.0000         0
%          0         0         0         0         0    1.0000
%          0   -7.3500    0.7875         0         0         0
%          0   73.5000  -33.0750         0         0         0
%          0   -6.0000    5.2143         0         0         0

% B =
% 
%          0
%          0
%          0
%     0.6071
%    -1.5000
%     0.4286

% syms x(t) theta1(t) theta2(t)
% syms x x_dot(t) theta1 theta1_dot(t) theta2 theta2_dot(t)
% 
% s = [x; 
%      diff(x); 
%      theta1; 
%      diff(theta1); 
%      theta2; 
%      diff(theta2)];
% 
% sdot = diff(s)
% sdot = diff(s) == A*s + B*u

%% ========= Control =========

% Eigenvalue Placement
% u = -kt -> x_dot = (A-B*k)*x
% poles : [0 0 -10.3827 10.3827 -4.0988 4.0988]
% eigs = [-73; -0.11+i; -0.11-i; -0.3; -6; -9];
% eigs = [-100; -0.1+i; -0.1-i; -0.3; -6; -9];
% K_eig_place = place(A,B,eigs);

% Linear Quadratic Regulator - LQR 
Q = diag([5 50 50 20 700 700]);
R = 1;       % K_lqr = [7.0711 -325.2388  -30.5972   17.7333  -50.7732  -82.2823]
% R = 0.1;   % K_lqr = [22.3607 -835.1590 -153.1012   55.0508 -141.9729 -248.9710]
% R = 0.001; % K_lqr = 1.0e+03 * [0.2236 -7.6870 -1.7568 0.5482 -1.3607 -2.4677]
K_lqr = lqr(A,B,Q,R);
% [K,S,P] = lqr(A,B,Q,R)
% A - State space matrix
% B - Input to state matrix
% Q - State cost weighted matrix
% R - Input cost weighted matrix
% K - Optimal gain
% S - Solution to the associated algebraic Riccatic Equation
% P - Poles of the closed loop system

% Open Loop Analysis

% Define the A and B matrices for the open-loop system (use the same matrices)
A_open = A;
B_open = B;
% Calculate the transfer function numerator and denominator for the open-loop system
[num_open, den_open] = ss2tf(A_open, B_open, C, D);

% Ensure that num_open and den_open are row vectors
num_open = num_open(1,:);
den_open = den_open(1,:);

% Create the transfer function for the open-loop system 
% 4 poles and 2 zeros
sys_open = tf(num_open, den_open);

% Closed loop Analysis

% Calculate closed-loop state-space matrices
A_closed = A - B*K_lqr; % This applies the control law (B*KK) to the unstable A matrix to stabilize it
B_closed = B;
C_closed = C;
D_closed = D;

% Calculate transfer function numerator and denominator for closed-loop system
[num_closed, den_closed] = ss2tf(A_closed, B_closed, C_closed, D_closed);

% Ensure that num_cl and den_cl are row vectors
num_closed = num_closed(1,:);
den_closed = den_closed(1,:);

% Create the transfer function for the closed-loop system
sys_closed = tf(num_closed, den_closed);

%% ========= Math =========

% ODE45 Notes

% ODE45 has nominally three inputs
% 1) [s_dot] = g(t,s)
% 2) tspan = [t0 tf] *wont impact accuracy
% 3) initial condition [x0, y0, v0, ..., ]

% [s] = [s1, s2, ... sn] - state variables of nth order for a differential equation

% Example: 
% x_triple_dot + x_double_dot + x_dot + x = 0
% [s] = [s1; = [x;
%        s2;    x_dot;
%        s3]    x_double_dot]

% x_triple_dot = - x_double_dot - x_dot - x 
% x_triple_dot = - s3 - s2 - s1 
% [s_dot] = [s1_dot; = [x_dot;        = [s2;           = [s2;
%            s2_dot;    x_double_dot;    s3;              s3;
%            s3_dot]    x_triple_dot]    x_triple_dot]   -s3-s2-s1]
 
% (m0 + m1 + m2)*x_double_dot + (m1*l1 + m2*L2)*cos(theta1)*theta1_double_dot + m2*l2*cos(theta2)*theta2_double_dot 
% - (m1*l1 + m2*L1)*sin(theta1)*(theta1_dot)^2 - m2*l2*sin(theta2)*(theta2_dot)^2 == u_t;

% (m1*(l1)^2 + m2*(L1)^2 + I1)*theta1_double_dot + (m1*l1 + m2*L1)*cos(theta1)*x_double_dot 
% + m2*L1*l2*cos(theta1-theta2)*theta2_double_dot + m2*L1*l2*sin(theta1-theta2)*(theta2_dot)^2 - g*(m1*l1 + m2*L1)*sin(theta1) == 0;

% m2*l2*cos(theta2)*x_double_dot + m2*L1*l2*cos(theta1-theta2)*theta1_double_dot + (m2*(l2)^2 + I2)*theta2_double_dot 
% - m2*L1*l2*sin(theta1-theta2)*(theta1_dot)^2 - m2*g*l2*sin(theta2) == 0;


% x0 = (0, theta1, theta2, 0, 0, 0)' - Initial Conditions
% x(t) = (0, 0, 0, 0, 0, 0)' - Pendulum up position

% syms x(t) theta1(t) theta2(t) u m0 m1 m2 L1 L2 l1 l2 I1 I2 g
% 
% eqn1_5 = (m0 + m1 + m2)*diff(x,2) + (m1*l1 + m2*L2)*cos(theta1)*diff(theta1,2) + m2*l2*cos(theta2)*diff(theta2,2) - (m1*l1 + m2*L1)*sin(theta1)*(diff(theta1,1))^2 - m2*l2*sin(theta2)*(diff(theta2,1))^2 == u;
% 
% eqn1_4 = (m1*(l1)^2 + m2*(L1)^2 + I1)*diff(theta1,2) + (m1*l1 + m2*L1)*cos(theta1)*diff(x,2) + m2*L1*l2*cos(theta1-theta2)*diff(theta2,2) + m2*L1*l2*sin(theta1-theta2)*(diff(theta2,1))^2 - g*(m1*l1 + m2*L1)*sin(theta1) == 0;
% 
% eqn1_3 = m2*l2*cos(theta2)*diff(x,2) + m2*L1*l2*cos(theta1-theta2)*diff(theta1,2) + (m2*(l2)^2 + I2)*diff(theta2,2) - m2*L1*l2*sin(theta1-theta2)*(diff(theta1,1))^2 - m2*g*l2*sin(theta2) == 0;
% 
% [sdot,s] = odeToVectorField(eqn1_4,eqn1_5,eqn1_3) % correct order

% M = matlabFunction(sdot,'vars',{'t','Y','u','m0','m1','m2','L1','L2','l1','l2','I1','I2','g'})

% x_dot             = diff(x);
% x_double_dot      = diff(x_dot);
% theta1_dot        = diff(theta1);
% theta1_double_dot = diff(theta1_dot);
% theta2_dot        = diff(theta2);
% theta2_double_dot = diff(theta2_dot);
% 
% eqn1_5 = (m0 + m1 + m2)*x_double_dot + (m1*l1 + m2*L2)*cos(theta1)*theta1_double_dot + m2*l2*cos(theta2)*theta2_double_dot - (m1*l1 + m2*L1)*sin(theta1)*(theta1_dot)^2 - m2*l2*sin(theta2)*(theta2_dot)^2 == u;
% 
% eqn1_4 = (m1*(l1)^2 + m2*(L1)^2 + I1)*theta1_double_dot + (m1*l1 + m2*L1)*cos(theta1)*x_double_dot + m2*L1*l2*cos(theta1-theta2)*theta2_double_dot + m2*L1*l2*sin(theta1-theta2)*(theta2_dot)^2 - g*(m1*l1 + m2*L1)*sin(theta1) == 0;
% 
% eqn1_3 = m2*l2*cos(theta2)*x_double_dot + m2*L1*l2*cos(theta1-theta2)*theta1_double_dot + (m2*(l2)^2 + I2)*theta2_double_dot - m2*L1*l2*sin(theta1-theta2)*(theta1_dot)^2 - m2*g*l2*sin(theta2) == 0;
% 
% [sdot,s] = odeToVectorField(eqn1_5,eqn1_4,eqn1_3)
% 
% M = matlabFunction(sdot,'vars',{'t','Y','u'})

%% ========= ODE Solver =========

% Initial Conditions - IC for simulation
% x          = -1; % Initial Horizontal postion of cart in [m]
% theta1_deg =  5; % Initial Angle between lower pendulum and the vertical in [deg]
% theta2_deg =  5; % Initial Angle between upper pendulum and the vertical in [deg]
% x_dot      =  0; % Initial cart velocity in [m/s]
% theta1_dot =  0; % Initial lower pendulum angular velocity in [m/s]
% theta2_dot =  0; % Initial upper pendulum angular velocity in [m/s]

x          =  0; % Initial Horizontal postion of cart in [m]
theta1_deg =  5; % Initial Angle between lower pendulum and the vertical in [deg]
theta2_deg =  5; % Initial Angle between upper pendulum and the vertical in [deg]
x_dot      =  0; % Initial cart velocity in [m/s]
theta1_dot =  0; % Initial lower pendulum angular velocity in [m/s]
theta2_dot =  0; % Initial upper pendulum angular velocity in [m/s]

xcom          = 0;
theta1com     = 0;
theta2com     = 0;
x_dotcom      = 0;
theta1_dotcom = 0;
theta2_dotcom = 0;

% Conversion from degrees to radians
theta1 = deg2rad(theta1_deg);
theta2 = deg2rad(theta2_deg);

t0           = 0;             % time initial to begin integration
tf           = 10;            % time final to finish integration
dt           = 0.01;          % time increment for integration
tspan        = t0:dt:tf; % time span for integration with time increment between initial and final time
state_values = [x theta1 theta2 x_dot theta1_dot theta2_dot]'; % Initial Conditions
y_com        = [xcom theta1com theta2com x_dotcom theta1_dotcom theta2_dotcom]'; % Commanded/desired state values to be achieved
u = @(y) -K_lqr*(y - y_com); % u - control input [scalar]

% u = @(y) -K_eig_place*(y - y_com); % u - control input [scalar]
[t, state_values] = ode45(@(t,y)dpic(y,g,L1,m0,m1,m2,l1,l2,I1,I2,u(y)), tspan, state_values);
% state_values % [2001x6] [rowxcol]

% Index into state_values and assign numerical values to state variables for plotting
%SV           = state_values;
x            = state_values(:,1);
theta1       = state_values(:,2);
theta2       = state_values(:,3);
x_dot        = state_values(:,4);
theta1_dot   = state_values(:,5);
theta2_dot   = state_values(:,6);

%% ========= Live simulation  =========

% theta1 = Y(2);
% theta2 = Y(3);
% theta1_dot = Y(5);
% theta2_dot = Y(6);

hf = figure(1);
movieVector = cell(1, ceil(numel(state_values)/8)); % Initialize cell array
% Calculate the maximum value of i
max_i = min(numel(state_values), 1000); % Limit to 1000 to avoid exceeding bounds
for i = 1:8:max_i
    % disp(['Current index: ', num2str(i)]);
    doubleIP_Animation(state_values(i,1), state_values(i,2), state_values(i,3), L1, L2);
    pause(0.001);
    movieVector{ceil(i/8)} = getframe(hf); % Store frame data in cell array
    hold off
end
%% =========  Save the movie ========= 
% myWriter = VideoWriter('double_IP.avi', 'Motion JPEG AVI');
% myWriter.Quality = 100;
% myWriter.FrameRate = 20;
% % Open the VideoWriter object, write the movie, and close the file
% open(myWriter);
% writeVideo(myWriter, cat(2, movieVector{:})); % Concatenate frames and write to video
% close(myWriter);

%% ========= Plot time domain results  =========

figure(2)

% Cart Position Plot
subplot(2,3,1)
plot(t,x,'LineWidth',1)
title('Cart Position')
xlabel('Time-t [s]')
ylabel('x [m]')
% xlim([0 10])
% ylim([0 3])
grid on

% Lower Pendulum Angle Plot
subplot(2,3,2)
plot(t,theta1,'LineWidth',1)
title('Lower Pendulum Angle')
xlabel('Time-t [s]')
ylabel('\theta_1 [rad]')
% xlim([0 10])
% ylim([-0.1 0.1])
grid on

% Upper Pendulum Angle Plot
subplot(2,3,3)
plot(t,theta2,'LineWidth',1)
title('Upper Pendulum Angle')
xlabel('Time-t [s]')
ylabel('\theta_2 [rad]')
% xlim([0 10])
% ylim([-0.1 0.1])
grid on

% Cart Velocity Plot
subplot(2,3,4)
plot(t,x_dot,'LineWidth',1)
title('Cart Velocity')
xlabel('Time-t [s]')
ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'latex')
% xlim([0 10])
% ylim([-1 1.5])
grid on

% Lower Pendulum Angular Velocity Plot
subplot(2,3,5)
plot(t,theta1_dot,'LineWidth',1)
title('Lower Pendulum Anglular Velocity')
xlabel('Time-t [s]')
ylabel('$\dot{\theta}_1$ [rad/s]', 'Interpreter', 'latex')
xlim([0 10])
% ylim([-0.15 0.1])
grid on

% Upper Pendulum Angular Velocity Plot
subplot(2,3,6)
plot(t,theta2_dot,'LineWidth',1)
title('Upper Pendulum Angular Velocity')
xlabel('Time-t [s]')
ylabel('$\dot{\theta}_2$ [rad/s]', 'Interpreter', 'latex')
xlim([0 10])
% ylim([-0.15 0.05])
grid on

x0     = 300;
y0     = 50;
width  = 900;
height = 700;
set(gcf,'position',[x0, y0, width, height])

% figure(3)
% plot(t,u)

%% Stability Analysis 

figure(4)
subplot(2,1,1);
rlocus(sys_open)
title('Open Loop Transfer Function Root Locus')
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
set(axIm,'String','Imaginary Component of S');
set(axRe,'String','Real Component of S');

subplot(2,1,2);
rlocus(sys_closed)
title('Closed Loop Transfer Function Root Locus')
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
set(axIm,'String','Imaginary Component of S');
set(axRe,'String','Real Component of S');

% Display the system's transfer functions

% Display the open-loop transfer function
disp('Open-Loop Transfer Function:');
sys_open  % Use a semi colon to suppress transfer function displaying in the console

% Display the transfer function for the closed-loop system
disp('Closed-Loop Transfer Function:');
sys_closed % Use a semi colon to suppress transfer function displaying in the console

% figure(5)
% U = zeros(size(t));
% for i = 1:length(state_values)
%     U(i) = u(state_values(i,:));
% end
% plot(t,U)
% title('Control Input');
% xlabel('Time-t [s]'); 
% ylabel('Control Input-u [N]');
% % xlim([0 1.45])
% ylim([-1.5 1.5])
% grid on;

%% Function definition for ODE45
function Y_dot = dpic(Y,g,L1,m0,m1,m2,l1,l2,I1,I2,u)

theta1     = Y(2);
theta2     = Y(3);
theta1_dot = Y(5);
theta2_dot = Y(6);

D = [m0+m1+m2                (m1*l1+m2*L1)*cos(theta1)     m2*l2*cos(theta2);
    (m1*l1+m2*L1)*cos(theta1) m1*(l1)^2+m2*(L1)^2+I1       m2*L1*l2*cos(theta1-theta2);
     m2*l2*cos(theta2)        m2*L1*l2*cos(theta1-theta2)  m2*(l2)^2+I2];

C = [0     -(m1*l1+m2*L1)*sin(theta1)*theta1_dot    -m2*l2*sin(theta2)*theta2_dot;
     0                      0                        m2*L1*l2*sin(theta1-theta2)*theta2_dot;
     0     -m2*L1*l2*sin(theta1-theta2)*theta1_dot   0];

G = [0;
   -(m1*l1+m2*L1)*g*sin(theta1);
    -m2*g*l2*sin(theta2)];

H = [1;
     0;
     0];

Z     = zeros(3); % Zeros matrix       [3x3] [rowxcol]
Z_col = [0 0 0]'; % Zero column vector [3x1] [rowxcol]
I     = eye(3);   % Identity matrix    [3x3] [rowxcol]
D_inv = inv(D);   % Inverse D matrix   [3x3] [rowxcol]

% D matrix 3x3 [rowxcol]
% C matrix 3x3 [rowxcol]
% G matrix 3x1 [rowxcol]
% H matrix 3x1 [rowxcol]

A = [Z  I;
     Z -D_inv*C]; % [3x3]*[3x3] [rowxcol]
B = [Z_col;
     D_inv*H];    % [3x3]*[3x1] [rowxcol]
L = [Z_col;
    -D_inv*G];    % [3x3]*[3x1] [rowxcol]

Y_dot = A*Y + B*u + L; % First order nonlinear system matrix model
end






