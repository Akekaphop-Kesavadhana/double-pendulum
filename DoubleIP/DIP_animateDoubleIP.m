function DIP_animateDoubleIP(t, z, p, varargin)
% DIP_animateDoubleIP  Smooth animation of a double-inverted pendulum on a cart,
%                      with a REPLAY button.
%
%   DIP_animateDoubleIP(t, z, p)                        % default options
%   DIP_animateDoubleIP(t, z, p, 'slowMo',1.5,'skip',2) % example
%
% INPUTS
%   t : 1xN   time vector
%   z : 6xN   [x; theta1; theta2; x_dot; theta1_dot; theta2_dot]
%   p : struct with fields  L1, L2   (link lengths, m)
%
% Name–Value options
%   'slowMo' (default 1)  :  playback factor (1 = real time, 2 = half speed, ...)
%   'skip'   (default 1)  :  plot every k-th sample

%---------------------------------------------------------------------
%  Parse options
%---------------------------------------------------------------------
opt.slowMo = 1.0;
opt.skip   = 1;
opt = parseNameValue(opt, varargin{:});

%---------------------------------------------------------------------
%  Pre-compute positions (for axis limits & fast access)
%---------------------------------------------------------------------
x   = z(1,:); theta1 = z(2,:); theta2 = z(3,:);
L1  = p.L1; L2 = p.L2;

p0  = [x; zeros(size(x))];                             % cart pivot
p1  = p0 + [ L1*sin(theta1);  L1*cos(theta1) ];        % 1st link tip  (0 rad ↑ )
p2  = p1 + [ L2*sin(theta2);  L2*cos(theta2) ];        % 2nd link tip

%---------------------------------------------------------------------
%  Figure and axes
%---------------------------------------------------------------------
hf = figure('Name','Double Inverted Pendulum – Swing-Up');
clf(hf); ax = axes('Parent',hf); hold(ax,'on'); axis(ax,'equal');
grid(ax,'on');
xlabel(ax,'X [m]'); ylabel(ax,'Y [m]');
title(ax,'Double Inverted Pendulum - Swing-Up');

pad = 0.3;
xlim(ax,[min([p0(1,:),p1(1,:),p2(1,:)])-pad,  max([p0(1,:),p1(1,:),p2(1,:)])+pad]);
ylim(ax,[min([p0(2,:),p1(2,:),p2(2,:)])-pad,  max([p0(2,:),p1(2,:),p2(2,:)])+pad]);
plot(ax,xlim(ax), [0 0],'k','LineWidth',2);          % ground

%---------------------------------------------------------------------
%  Graphics objects (handles) - will be updated each frame
%---------------------------------------------------------------------
cartW = 0.3; cartH = 0.15;
cart   = rectangle(ax,'Position',[x(1)-cartW/2, -cartH/2, cartW, cartH], ...
                   'Curvature',0.1,'FaceColor',[0.8 0 0],'EdgeColor','none');
link1  = plot(ax,[p0(1,1) p1(1,1)],[p0(2,1) p1(2,1)],'b','LineWidth',3);
link2  = plot(ax,[p1(1,1) p2(1,1)],[p1(2,1) p2(2,1)],'b','LineWidth',3);
joint1 = plot(ax,p1(1,1),p1(2,1),'ro','MarkerFaceColor','r','MarkerSize',6);
joint2 = plot(ax,p2(1,1),p2(2,1),'ro','MarkerFaceColor','r','MarkerSize',6);

%---------------------------------------------------------------------
%  REPLAY button
%---------------------------------------------------------------------
uicontrol('Parent',hf, 'Style','pushbutton', 'String','Replay', ...
          'Units','normalized', 'Position',[0.88 0.02 0.10 0.06], ...
          'FontWeight','bold', 'Callback',@(~,~) runAnimation);

%---------------------------------------------------------------------
%  Run once automatically, then every time user presses Replay
%---------------------------------------------------------------------
runAnimation(); % initial play-through

% ======================== nested functions ==========================
    function runAnimation()
        skip   = max(1, round(opt.skip));
        slowMo = max(0.01,opt.slowMo);
        N      = length(t);

        for k = 1:skip:N
            % --- update cart
            cart.Position(1) = x(k) - cartW/2;

            % --- update links
            link1.XData = [p0(1,k) p1(1,k)];  link1.YData = [p0(2,k) p1(2,k)];
            link2.XData = [p1(1,k) p2(1,k)];  link2.YData = [p1(2,k) p2(2,k)];

            % --- update joints
            joint1.XData = p1(1,k); joint1.YData = p1(2,k);
            joint2.XData = p2(1,k); joint2.YData = p2(2,k);

            drawnow limitrate

            % pacing
            if k>skip
                dt = t(k) - t(k-skip);
                pause(slowMo*dt);
            end
        end
    end
end
% ====================== helper: parse NV pairs =====================
function opt = parseNameValue(opt, varargin)
for ii = 1:2:length(varargin)
    name  = varargin{ii};
    value = varargin{ii+1};
    if isfield(opt,name)
        opt.(name) = value;
    else
        error('Unknown option "%s"', name);
    end
end
end
