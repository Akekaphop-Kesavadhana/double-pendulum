function u_app = DIP_motor_system_dynamics(velocity,u_com)

% INPUTS:
%     velocity = x_dot 
%     u_com = The commanded force to be applied to the cart   
% OUTPUTS: 
%     u_applied = The achieved force on the cart

wm   = 4;               % Motor angular velocity   [rad/s]
wp   = 2;               % Pulley angular velocity  [rad/s]
dp   = 0.0254;          % Pulley diameter          [m]
R    = 1.9;             % Motor winding resistance [Ohm]
Ktau = 0.104511;        % Motor torque Constant    [N*m/A]
Ke   = 0.105042;        % Motor back-EMF Constant  [V/(rad/s)]
% Ir = 4.801855e-5;     % Rotor inertia Constant   [N*m/s^2] NOT USED
G    = wm/wp;           % Gear ratio 
                        % G>1 pulley spins slower than motor 
                        % G<1 pulley spins faster than motor

% dp   = 0.01;    % pulley radius [m]
% R    = 1.0;     % winding resistance [Ω]
% Ktau = 0.163;   % torque constant [N·m/A]
% Ke   = 0.163;   % back-EMF constant [V·s/rad]
% G    = 1;       % gear ratio
vBus = 12;        % Supply limit [V]

Kff  = (R*dp)/(G*Ktau);
Kemf = (G*Ke)/dp;
KVF  = (G*Ktau)/(dp*R);
Kvf  = (G^2*Ktau*Ke)/(dp^2*R);

% --- one-line scalar computation ---
% voltage   = max( min(Kff*u_com + Kemf*velocity,  vBus), -vBus ); min and max functions cause solver instability
voltage   = Kff*u_com + Kemf*velocity;
u_app = KVF*voltage - Kvf*velocity;
end