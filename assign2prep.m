%% Get data.
% This script does some preprocessing of the experimental data recorded
% from the single direct drive motor.

%% Closed-loop data

load('assign2cl.mat');

% Note that the variables i this file are:
%  t = Nx1 vector with the time of the N samples.
%  y = Nx3 matrix with for each sample:
%          y(:,1) -> measured angle (n [rad]
%          y(:,2) -> applied torque in [Nm]
%          y(:,3) -> reference angle in [rad]

%% Number of samples and sample time

N = length(t);
t_s = (t(end)-t(1))/(N-1);

%% Extract data
u = y(:,2);     % Torque should be in 2nd column
r = y(:,3);     % Reference angle in 3rd column
y = y(:,1);     % Measured angle in 1st column. Note that y is overwritten!!!

%% Filter torque with moving average
frange = 30;    % Width of the moving average filter in samples (minus 1)
uf = zeros(size(u));
for ind = 1:length(u)
    if ind<frange/2+1
        uf(ind) = mean(u(1:ind));   % Less samples at start of data
    elseif ind>length(u)-frange/2
        uf(ind) = mean(u(ind:end)); % Less sampkes at end of data 
    else
        uf(ind) = mean(u(ind-frange/2:ind+frange/2));
    end
end

%% Compute angular velocity and direction
% Some low-pass filtering is applied in the computation of the velocity by
% evaluating the difference in a longer time interval of +/- dfilt*t_s:
dfilt = 7;
yv = [zeros(dfilt,1); y(2*dfilt+1:end)-y(1:end-2*dfilt); zeros(dfilt,1)] ...
     /(2*dfilt*t_s);    % [ rad/s]
yd = sign(yv);

%% Set range for data to be used in fit. 
% These are the reference angles >=0 and <2 pi

range = find(r>=0 & r<2*pi);

t = t(range);
r = r(range);
u = u(range);
y = y(range);
uf= uf(range);
yd= yd(range);
yv= yv(range);
