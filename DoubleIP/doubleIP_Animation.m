function animation =  doubleIP_Animation(x,th,th2,len1,len2) 

% Parameter x - cart position, th - theta1, th2 - theta2

persistent j
persistent pjx
persistent pjy
    if isempty(j)
        j = 0;
    end
    j = j+1;
    
W   = 0.2;  % width of cart
H   = 0.08; % height of Cart
wr  = 0.04; % right wheel 
L1  = len1; % length of first pendulum 
L2  = len2; % length of second pendulum 

% Position coordination of the first pendulum
y   = H/2+wr/2;
w1x = x -0.9*W/2;
w1y = 0;
w2x = x+0.9*W/2-wr;
w2y = 0;

% Position coordination of the second pendulum
y2   = H/2+wr/2;
w1x2 = x -0.9*W/2;
w1y2 = 0;
w2x2 = x+0.9*W/2-wr;
w2y2 = 0;


% Position of 1st pendulum 
px = x - L1*sin(th);
py = y + L1*cos(th);
% Position of 2nd pendulum 
pxx = px - L2*sin(th2);
pyy = py + L2*cos(th2);

pjx(j) = pxx;
pjy(j) = pyy;

base        = plot([-4 4],[0 0],'k','LineWidth',2); % Base line/ground
hold on;
cart        = rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',0.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1]);
left_wheel  = rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
right_wheel = rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
    
pendulum    = plot([x px],[y py],'b','LineWidth',2.5);                        % First Pendulum rod 
pendulum2   = plot([px pxx],[py pyy],'b','LineWidth',2.5);                    % Second Pendulum rod 
p_cir       = viscircles([px py],0.02,'Color',[1 0.1 0.1],'LineWidth',2.5);   % Pendulum Trajectory(Circle one)
p_cir2      = viscircles([pxx pyy],0.02,'Color',[1 0.1 0.1],'LineWidth',2.5); % Pendulum Trajectory(Circle two)
p_cir1      = viscircles([x y],0.02,'Color','w','LineWidth',0.2);             % center of Cart

% line_traj = plot(pjx(1:j),pjy(1:j), 'm--','LineWidth',1);  % Pendulum Trajectory (Line)
    xlabel('X [m]');
    ylabel('Y [m]');
    title('Double Inverted Pendulum: SwingUp Control')
    axis(gca,'equal');
    xlim([-4 4]);
    ylim([-1.5 1.5]);
    x0     = 300;
    y0     = 50;
    width  = 900;
    height = 700;
    set(gcf,'position',[x0, y0, width, height])
    grid on;
    % drawnow
   