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

% Citation :
% Crowe-Wright, Ian J P. "Control Theory: The Double Pendulum Inverted on a Cart." (2018). https://digitalrepository.unm.edu/math_etds/132

% Constant initializations and definitions
g  = 9.8;               % Gravitational constant in [m/s^2]
m0 = 1;                 % Mass of cart in [kg]
m1 = 0.5;                 % Mass of lower pendulum link in [kg]
m2 = 0.5;                 % Mass of upper pendulum link in [kg]
L1 = 0.5;               % Length of lower pendulum link in [m]
L2 = 0.5;               % Length of upper pendulum link in [m]
l1 = (1/2)*L1;          % Distance between the base of first pend. link and its center of mass in [m]
l2 = (1/2)*L2;          % Distance between the base of second pend. link and its center of mass in [m]
I1 = (1/12)*m1*(L1)^2;  % Moment of inertia of first pendulum link in [kg*(m)^2]
I2 = (1/12)*m2*(L2)^2;  % Moment of inertia of second pendulum link in [kg*(m)^2]
n  = 2;
con = [I1,I2,m0,m1,m2,L1,L2,l1,l2,n,g]';

%% ========= Linearized state space form and matrix constituents =========
% This linear model is only valid for the condition when both pendulums are upright and
% possess angles of approximately +-39 degrees relative to the vertical 

p = 1/(4*m0*m1 + 3*m0*m2 + (m1)^2 + m1*m2);

a42 = -(3/2)*p*(2*(m1)^2 + 5*m1*m2 + 2*(m2)^2)*g;
a43 =  (3/2)*p*m1*m2*g;
a52 =  (3/2)*(p/L1)*(4*m0*m1 + 8*m0*m2 + 4*(m1)^2 + 9*m1*m2 + 2*(m2)^2)*g;
a53 = -(9/2)*(p/L1)*(2*m0*m2 + m1*m2)*g;
a62 = -(9*p*g*(2*m0+m1)*(m1+2*m2))/(2*L2);  
a63 =  (3/2)*((p*g)/L2)*((m1)^2 + 4*m0*m1 + 12*m0*m2 + 4*m1*m2); 

b4 = p*(4*m1 + 3*m2);
b5 =  -((3*p)/L1)*(2*m1 + m2);
b6 =   (3*p*m1)/L2; % m2 changed to m1 in b6

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

%% ========= Linearized state space form and matrix constituents (Lower pend stabilization once swung up) =========

% a1 = m1*l1+L1*(m2+n2);
% a2 = -(c1+c2)/(I1+m2*(L1)^2+n2*(L1)^2);
% 
% b1 = -(m1*l1+m2*L1+n2*L1)/(I1+m2*(L1)^2+n2*(L1)^2);
% 
% A1 = [0  1  0 0;
%       a1 a2 0 0;
%       0  0  0 1;
%       0  0  0 0];
% 
% B1 = [0
%       b1
%       0
%       1];

%% ========= Linearized state space form and matrix constituents (Lower + Upper pend stabilization once upper pend is swung up) =========

% Itheta1 = I1 + m2*(L1)^2 + n2*(L1)^2;
% Itheta2 = I2 + Jn1 + Jn2;
% M1      = m2*l2*L1;
% M2      = m1*l1 + m2*L1 + n2*L1;
% M3      = m2*l2;
% 
% a221 = -(M2*Itheta2*g)/(Itheta2*Itheta1-(M1)^2);
% a222 = -(M1*c2+Itheta2*(c1+c2))/(Itheta2*Itheta1-(M1)^2);
% a223 =  (M1*M3*g)/(Itheta2*Itheta1-(M1)^2);
% a224 =  (Itheta2*c1 + M1*c2)/(Itheta2*Itheta1-(M1)^2);
% a241 =  (M1*M2*g)/(Itheta2*Itheta1-(M1)^2);
% a242 =  (Itheta1*c2+M1*(c1+c2))/(Itheta2*Itheta1-(M1)^2);
% a243 = -(M3*Itheta1*g)/(Itheta2*Itheta1-(M1)^2);
% a244 = -(Itheta1*c2+M1*c1)/(Itheta2*Itheta1-(M1)^2);
% 
% b21 = (M1*M3-M2*Itheta2)/(Itheta2*Itheta1-(M1)^2);
% b41 = (M1*M2-M3*Itheta1)/(Itheta2*Itheta1-(M1)^2);
% 
% 
% A2 = [0    1    0    0    0 0;
%       a221 a222 a223 a224 0 0;
%       0    0    0    1    0 0;
%       a241 a242 a243 a244 0 0;
%       0    0    0    0    0 1;
%       0    0    0    0    0 0];
% 
% B2 = [0;
%       b21;
%       0;
%       b41;
%       0;
%       1];


%% ========= Control =========

% Linear Quadratic Regulator - LQR 
Q = diag([50 50 50 20 700 700]);
R = 1;       
K_lqr = lqr(A,B,Q,R);

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

%% ========= ODE Solver =========

% Initial Conditions - IC for simulation
x          = 0;  % Initial Horizontal postion of cart in [m]
theta1_deg = 180;   % Initial Angle between lower pendulum and the vertical in [deg]
theta2_deg = 180;   % Initial Angle between upper pendulum and the vertical in [deg]
x_dot      = 0;   % Initial cart velocity in [m/s]
theta1_dot = 0;   % Initial lower pendulum angular velocity in [m/s]
theta2_dot = 0;   % Initial upper pendulum angular velocity in [m/s]

% Desired states
x_com          = 0;   % Desired Horizontal postion of cart in [m]
theta1_deg_com = 0;   % Desired Angle between lower pendulum and the vertical in [deg]
theta2_deg_com = 0;   % Desired Angle between upper pendulum and the vertical in [deg]
x_dot_com      = 0;   % Desired cart velocity in [m/s]
theta1_dot_com = 0;   % Desired lower pendulum angular velocity in [m/s]
theta2_dot_com = 0;   % Desired upper pendulum angular velocity in [m/s]

% Conversion from degrees to radians
theta1     = deg2rad(theta1_deg);
theta2     = deg2rad(theta2_deg);
theta1_com = deg2rad(theta1_deg_com);
theta2_com = deg2rad(theta2_deg_com);

t0 = 0;             % time initial to begin integration
tf = 30;            % time final to finish integration
tspan = t0:0.01:tf; % time span for integration with time increment between initial and final time

IC    = [x     theta1     theta2     x_dot     theta1_dot     theta2_dot]'; % Initial Conditions
y_com = [x_com theta1_com theta2_com x_dot_com theta1_dot_com theta2_dot_com]'; % Commanded/desired state values to be achieved

% u = @(y) -K_lqr*(y - y_com); % u - control input [scalar]
u = @(y) controller(y,con,y_com);
[t, y] = ode45(@(t,y)dpic(y,g,L1,m0,m1,m2,l1,l2,I1,I2,u(y)), tspan, IC);

% Index into state_values and assign numerical values to state variables for plotting
x            = y(:,1);
theta1       = y(:,2);
theta2       = y(:,3);
x_dot        = y(:,4);
theta1_dot   = y(:,5);
theta2_dot   = y(:,6);

%% ========= Plot time domain results  =========

% Visual simulation 
hf = figure(1);
movieVector = cell(1, ceil(numel(y)/8)); % Initialize cell array

% Calculate the maximum value of i
max_i = min(numel(y), 1000); % Limit to 1000 to avoid exceeding bounds

for i = 1:8:max_i
    % disp(['Current index: ', num2str(i)]);
    doubleIP_Animation(y(i,1), y(i,2), y(i,3), L1, L2);
    pause(0.001);
    movieVector{ceil(i/8)} = getframe(hf); % Store frame data in cell array
    hold off
end

%% Save the movie
% myWriter = VideoWriter('double_IP.avi', 'Motion JPEG AVI');
% myWriter.Quality = 100;
% myWriter.FrameRate = 20;
% 
% % Open the VideoWriter object, write the movie, and close the file
% open(myWriter);
% writeVideo(myWriter, cat(2, movieVector{:})); % Concatenate frames and write to video
% close(myWriter);

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

%% Stability Analysis 

figure(3)
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

A = [Z  I;
     Z -D_inv*C];  % [3x3]*[3x3] [rowxcol]
B = [Z_col;
     D_inv*H];     % [3x3]*[3x1] [rowxcol]
L = [Z_col;
     -D_inv*G];    % [3x3]*[3x1] [rowxcol]

Y_dot = A*Y + B*u + L; % First order nonlinear system matrix model
end

function utot = controller(y,CON,yr)
    x            = y(1);
    theta1       = y(2);
    theta2       = y(3);
    x_dot        = y(4);
    theta1_dot   = y(5);
    theta2_dot   = y(6);

    xr           = yr(1);
    theta1r      = yr(2);
    theta2r      = yr(3);
    x_dotr       = yr(4);
    theta1_dotr  = yr(5);
    theta2_dotr  = yr(6);
    
    I1 = CON(1);
    I2 = CON(2);
    m0 = CON(3);
    m1 = CON(4);
    m2 = CON(5);
    L1 = CON(6);
    L2 = CON(7);
    l1 = CON(8);
    l2 = CON(9);
    n  = CON(10);
    g  = CON(11);

%% ========= Lower pend swing up + Stabilization =========

    A1 = [0 0 1 0;
         0 0 0 1;
         0 -1.6259  0 0;
         0  63.4358 0 0];

    B1 = [0 0 0.6454 -2.2979]';

    y1  = [x theta1 x_dot theta1_dot]';
    y1r = [xr theta1r x_dotr theta1_dotr]';

    % Linear Quadratic Regulator - LQR 
    Q1 = diag([200 700 200 700]);
    R1 = 0.01;       
    K_lqr1  = lqr(A1,B1,Q1,R1);
    u1_lqr = -K_lqr1*(y1 - y1r);

    % Swing-up Controller based on energy
    E1 = (1/2)*m1*(l1)^2*theta1_dot^2 + m1*g*l1*(cos(theta1)-1);
    E0 = 0;  % Energy threshold (can be adjusted)
    u1_swing = sat(n*g*siggn((E1-E0)*(theta1_dot*cos(theta1))),-7,7);

    % Control logic (choose based on condition)
    if abs(theta1) < deg2rad(10)  % If close to upright position, use LQR
        u_control = u1_lqr;
    else  % Else, use swing-up controller
        u_control = u1_swing;
    end

%% ========= stabilization of lower pend while simulataneously swinging up the upper pend =========
    % E2 = (1/2)*I2*theta2_dot + m2*g*l2*(cos(theta2)-1);

    p = 1/(4*m0*m1 + 3*m0*m2 + (m1)^2 + m1*m2);

    a42 = -(3/2)*p*(2*(m1)^2 + 5*m1*m2 + 2*(m2)^2)*g;
    a43 =  (3/2)*p*m1*m2*g;
    a52 =  (3/2)*(p/L1)*(4*m0*m1 + 8*m0*m2 + 4*(m1)^2 + 9*m1*m2 + 2*(m2)^2)*g;
    a53 = -(9/2)*(p/L1)*(2*m0*m2 + m1*m2)*g;
    a62 = -(9*p*g*(2*m0+m1)*(m1+2*m2))/(2*L2);  
    a63 =  (3/2)*((p*g)/L2)*((m1)^2 + 4*m0*m1 + 12*m0*m2 + 4*m1*m2); 

    b4 = p*(4*m1 + 3*m2);
    b5 =  -((3*p)/L1)*(2*m1 + m2);
    b6 =   (3*p*m1)/L2; % m2 changed to m1 in b6
    
   % x    x_dot  theta1 theta1_dot theta2 theta2_dot
    A2 = [0    0      0      1     0      0;        % System matrix       s = [x
          0    0      0      0     1      0;                              %    x_dot
          0    0      0      0     0      1;                              %    theta1
          0    a42    a43    0     0      0;                              %    theta1_dot
          0    a52    a53    0     0      0;                              %    theta2
          0    a62    a63    0     0      0];                             %    theta2_dot]

    B2 = [0;
          0;
          0;
          b4;
          b5;
          b6];
    y2  = [x theta1 theta2 x_dot theta1_dot theta2_dot]';
    yr2 = [xr theta1r theta2r x_dotr theta1_dotr theta2_dotr]';
    
    Q2 = diag([50 50 50 20 700 700]);
    R2 = 1;       
    K_lqr2  = lqr(A2,B2,Q2,R2);
    u2_lqr = -K_lqr2*(y2 - yr2);
    % u_control = u2_lqr;
    % 
    % % Control logic (choose based on condition)
    % if abs(theta1) <= (deg2rad(10)) && abs(theta2) <= (deg2rad(10))  % If close to upright position, use LQR
    %     u_control = u2_lqr;
    % else  % Else, use swing-up controller
    %     u_control = 1;%u2_swing;
    % end

    % if (theta2_dot.*sin(theta2) < 0) & (cos(theta2) < cos(pi/6))
    %     uz2 = ua2.*sign(E2.*theta2_dot.*cos(theta2_dot));
    % elseif theta2_dot.*sin(theta2) >= 0 & cos(theta2) < cos(pi/6)
    %     uz2 = 0;
    % else, cos(theta2) >= cos(pi/6);
    %     uz2 = 0;
    % end
    % 
    % A11 =  Itheta1;
    % A22 =  Itheta2;
    % A12 =  M1.*cos(theta1-theta2);
    % f11 =  M1.*(theta2_dot).^2.*sin(theta1-theta2)+(c1+c2).*theta1_dot-c2.*theta2_dot-M2.*g.*sin(theta1);
    % f12 =  M2.*cos(theta1);
    % f21 = -M1.*(theta1_dot).^2.*sin(theta1-theta2)-c2.*(theta1_dot-theta2_dot)-M3.*g.*sin(theta2_dot);
    % f22 =  M3.*cos(theta2);
    % 
    % A33 = A11.*A22-A12.*A12 + (A12.*f22 - A22.*f12).*L1.*cos(theta1);
    % u2  = (uz2 + L1.*(theta1_dot).^2.*sin(theta1).*(A11.*A22 - A12.*A12) - ((A12.*f21 - A22.*f11).*L1.*cos(theta1)))/A33;

     % Apply boundary conditions on the cart's position and velocity
    max_pos =  1.5;  % Maximum cart position (m)
    min_pos = -1.5;  % Minimum cart position (m)
    k       =  2;

    % Ensure the cart's position stays within bounds by adjusting the control input
    if x >= max_pos  % If the position exceeds max
        u_control = min(u_control, -k*x_dot);  % Apply no control to move back to limit
    elseif x <= min_pos  % If the position exceeds min
        u_control = max(u_control, k*x_dot);  % Apply no control to move back to limit
    end

    utot = u_control;
end

function output = siggn(input)
    if input >= 0
        output = 1;
    else
        output = -1;
    end
end

function output = sat(input, lower_limit, upper_limit)
    output = min(max(input, lower_limit), upper_limit);
end




