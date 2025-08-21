close all; clear; clc; 

addpath ../../

L1 = 0.75;
L2 = 0.5;
m1 = 100000;
m2 = 0.5;

p.m0 = 2.0;               % Mass of cart [kg] [data type: struct array]
p.m1 = m1;                % Mass of lower pendulum link [kg] 
p.m2 = m2;                % Mass of upper pendulum link [kg]
p.L1 = L1;                % Length of lower pendulum link [m]
p.L2 = L2;                % Length of upper pendulum link [m]
p.l1 = (1/2)*L1;          % Distance between the base of first pend. link and its center of mass [m]
p.l2 = (1/2)*L2;          % Distance between the base of second pend. link and its center of mass [m]
p.I1 = (1/12)*m1*(L1)^2;  % Moment of inertia of first pendulum link [kg*(m)^2]
p.I2 = (1/12)*m2*(L2)^2;  % Moment of inertia of second pendulum link [kg*(m)^2]
p.g  = 9.81;              % Gravitational constant [m/s^2]

% dist     = 2;    % How far must the cart translate during its swing-up [m]
% maxForce = 80;   % Maximum actuator force [N]

% dist     = 0.5;  % How far must the cart translate during its swing-up [m]
maxForce   = 88;   % Maximum actuator force [N]
% Vmax     = 24;   % Maximum voltage applied to the motor [V]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)(DIP_dynamics_TrajOpt(x,u,p)); % [1x1 struct] System equations of motion
% problem.func.pathObj = @(t,x,u)(ones(size(t)));              % [1x1 struct] Swing-up Time cost function
% problem.func.pathObj = @(t,x,u)(u.^2);                       % [1x1 struct] Force-squared cost function
% problem.func.pathCst = @(t,x,u) voltagePath(t,x,u,Vmax);     % [1x1 struct]

% Weights
w_x = 30;    % Penalize cart displacement
w_u = 30;    % Penalize control effort 
xf  = 0;     % Final x (cart) position [m]

% path objective: ∫ [ w_x*(x–xi).^2  +  w_u*u.^2 ] dt
problem.func.pathObj = @(t,x,u)(w_x*(x(1,:) - xf).^2 + w_u*(u).^2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up problem bounds                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low   = 0.01;
problem.bounds.finalTime.upp   = inf;

% ------------------ Initial boundary values ------------------ %
xi           = 0;   % Initial x (cart) position [m]
theta1_deg_i = 180; % Initial theta1 angle [deg]
theta2_deg_i = 180; % Initial theta2 angle [deg]
theta1_rad_i = deg2rad(theta1_deg_i); theta2_rad_i = deg2rad(theta2_deg_i);
problem.bounds.initialState.low = [xi; theta1_rad_i; theta2_rad_i; 0; 0; 0];
problem.bounds.initialState.upp = [xi; theta1_rad_i; theta2_rad_i; 0; 0; 0];

% ------------------ Final boundary values ------------------ %
% xf         = 0; % Final x (cart) position [m]
theta1_deg_f = 0; % Final theta1 angle [deg]
theta2_deg_f = 0; % Final theta2 angle [deg]
theta1_rad_f = deg2rad(theta1_deg_f); theta2_rad_f = deg2rad(theta2_deg_f);
problem.bounds.finalState.low = [xf;theta1_rad_f;theta2_rad_f;0;0;0];
problem.bounds.finalState.upp = [xf;theta1_rad_f;theta2_rad_f;0;0;0];

% ------------------ Path boundaries ------------------ %
MaxU                     = 0.7; % +-Maximum travel distance [m]
problem.bounds.state.low = [-MaxU;-2*pi;-2*pi;-inf;-inf;-inf];
problem.bounds.state.upp = [MaxU; 2*pi; 2*pi; inf; inf; inf];

% ------------------ Maximum commanded force ------------------ %
problem.bounds.control.low = -maxForce;
problem.bounds.control.upp = maxForce;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time    = [0,2];
problem.guess.state   = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [0,0];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options(1).nlpOpt = optimset(...
    'Display','iter',...
    'TolFun',1e-8,...
    'MaxFunEvals',2e6);
problem.options(1).method = 'trapezoid';
problem.options(1).trapezoid.nGrid = 10;

problem.options(2).nlpOpt = optimset(...
    'Display','iter',...
    'TolFun',1e-8,...
    'MaxFunEvals',2e6);
problem.options(2).method = 'trapezoid';
problem.options(2).trapezoid.nGrid = 30;

% Voltage Path Constraint
% problem.options(3).nlpOpt = optimset(...
%     'Display','iter',...
%     'TolFun',1e-8,...
%     'MaxFunEvals',2e6);
% problem.options(2).method = 'trapezoid';
% problem.options(3).trapezoid.nGrid = 60;   

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Post Processing and Plots                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% ------------------ Unpack the simulation ------------------ %
t = linspace(soln(end).grid.time(1), soln(end).grid.time(end), 150);
z = soln(end).interp.state(t);
u = soln(end).interp.control(t);

% ------------------ Plot the state and control against time % ------------------ %
figure(1); clf;
DIP_plot_TrajOpt(t,z,u,p);

% ------------------ Visualise the piece-wise polynomial that optimTraj uses internally  ------------------ %
DIP_plotInterpolatedTraj(soln); % default 1000-point grid

% ------------------ Draw Trajectory ------------------ %
% [p1,p2] = cartPoleKinematics(z,p);
% figure(2); clf;
% nFrame = 9;  %Number of frames to draw
% drawCartPoleTraj(t,p1,p2,nFrame);

% ------------------ Display error in the collocation constraint between grid points ------------------ %
% if strcmp(soln(end).problem.options.method,'trapezoid') || strcmp(soln(end).problem.options.method,'hermiteSimpson')
%     % Then we can plot an estimate of the error along the trajectory
%     figure(5); clf;
% 
%     % NOTE: the following commands have only been implemented for the direct
%     % collocation(trapezoid, hermiteSimpson) methods, and will not work for
%     % chebyshev or rungeKutta methods.
%     cc = soln(end).interp.collCst(t);
% 
%     subplot(2,2,1);
%     plot(t,cc(1,:))
%     title('Collocation Error:   dx/dt - f(t,x,u)')
%     ylabel('d/dt cart position')
% 
%     subplot(2,2,3);
%     plot(t,cc(2,:))
%     xlabel('time')
%     ylabel('d/dt pole angle')
% 
%     idx = 1:length(soln(end).info.error);
%     subplot(2,2,2); hold on;
%     plot(idx,soln(end).info.error(1,:),'ko');
%     title('State Error')
%     ylabel('cart position')
% 
%     subplot(2,2,4); hold on;
%     plot(idx,soln(end).info.error(2,:),'ko');
%     xlabel('segment index')
%     ylabel('pole angle');
% end

% ------------------ Animation ------------------ %
t  = linspace(soln(end).grid.time(1), soln(end).grid.time(end), 300);
z  = soln(end).interp.state(t);            % 6xN 
DIP_animateDoubleIP(t, z, p, 'slowMo', 1); % 1.2 = 20 % slower than real time

% function [c, ceq] = voltagePath(~, x, u, Vmax)
% % voltagePath  Path inequality that limits motor voltage to +-Vmax.
% %
% %   [c, ceq] = voltagePath(~, x, u, Vmax)
% %
% %   c   : (2xN) inequality rows must satisfy  c <= 0
% %   ceq : [] (no equalities)
% 
%     x_dot = x(4,:);                              % cart velocity row-vector
%     v     = DIP_motor_system_dynamics(x_dot, u); % raw motor voltage [V]
% 
%     c   = [ v - Vmax ;          %  v ≤ +Vmax
%            -v - Vmax ];         % -v ≤ +Vmax
%     ceq = [];                   % no equality constraints
% end