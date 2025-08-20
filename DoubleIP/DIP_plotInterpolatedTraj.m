function DIP_plotInterpolatedTraj(soln, nSample)
% DIP_plotInterpolatedTraj  Plot all states and the control from an
%                           optimTraj solution, showing interpolant vs. nodes.
%
%   DIP_plotInterpolatedTraj(soln)             % default 1000 samples
%   DIP_plotInterpolatedTraj(soln, nSample)    % specify dense-grid length
%
% INPUT
%   soln     : structure returned by optimTraj (use last element if multi-pass)
%   nSample  : (optional) number of points for dense sampling (default 1000)

% -------------------------------------------------------------------------
%  Defaults
% -------------------------------------------------------------------------
if nargin < 2 || isempty(nSample)
    nSample = 1000;
end

% grab the last (finest) solution pass if the user passed the whole array
if numel(soln) > 1
    soln = soln(end);
end

% -------------------------------------------------------------------------
%  Dense sampling of the interpolant
% -------------------------------------------------------------------------
tDense = linspace(soln.grid.time(1), soln.grid.time(end), nSample);

zDense = soln.interp.state(tDense);    % [nState x nSample]
uDense = soln.interp.control(tDense);  % [nCtrl  x nSample]

nState = size(zDense,1);
hasControl = ~isempty(soln.grid.control);

% -------------------------------------------------------------------------
%  Prepare figure
% -------------------------------------------------------------------------
nRow = nState + double(hasControl);
stateNames = compose('state %d', 1:nState);
ctrlLabel  = 'control  u';

figure('Name','OptimTraj - Interpolated Trajectories'); clf;

% -------------------------------------------------------------------------
%  Plot each state
% -------------------------------------------------------------------------
for i = 1:nState
    subplot(nRow,1,i); hold on; grid on;
    plot(tDense, zDense(i,:), 'b-', 'LineWidth',1.4);
    plot(soln.grid.time, soln.grid.state(i,:), 'ro', 'MarkerSize',4);
    axis tight
    ylabel(stateNames{i});
    if i == 1
        title('Interpolated trajectory (blue) vs. collocation nodes (red)');
    end
end

% -------------------------------------------------------------------------
%  Plot control (if present)
% -------------------------------------------------------------------------
if hasControl
    subplot(nRow,1,nRow); hold on; grid on;
    plot(tDense, uDense, 'b-', 'LineWidth',1.4);
    plot(soln.grid.time, soln.grid.control, 'ro', 'MarkerSize',4);
    axis tight
    ylabel(ctrlLabel);
    xlabel('time  [s]');
else
    xlabel(subplot(nRow,1,nRow-1),'time  [s]');
end
end
