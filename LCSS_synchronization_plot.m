% Synchronization of a dynamical network of Lorenz oscillators
% ------------------------------------------------------------
% This script simulates the synchronization of 10 Lorenz oscillators 
% over two time intervals with different connectivity digraphs (G1 and G2).
% Synchronization error data is saved in "sync_data.xlsx".

clc; clear; close all;

%% Lorenz oscillator parameters
sigma = 10;
rho   = 25;
beta  = 8/3;

% Coupling strength (analytically derived)
a = -sigma + (beta*(beta+1)*(rho+sigma)^2) / (16*(beta-1));

%% Define connectivity digraphs
% Connectivity for first time interval
tail1 = [1 2 3 3 4 8 8 7 8 6 7 9 10 5 5 7];
head1 = [2 3 1 4 1 1 3 3 6 7 8 10 5 9 8 9];
G1    = digraph(tail1, head1);

% Connectivity for second time interval
tail2 = [1 2 2 3 3 4 1 4 2 3 4 5 5 6 7 8 8 8 7 9 10 1];
head2 = [2 1 3 1 4 1 5 5 6 6 7 6 7 8 5 6 7 9 9 10 9 10];
G2    = digraph(tail2, head2);

N         = 10;  % Number of oscillators (nodes)
numStates = 3;   % State dimension of each oscillator

%% Simulation settings
data_length1 = 25;
data_length2 = 25;
t_end1 = 1.5;                  % End of first interval
t_end2 = 1.5;                  % End of second interval
tspan1 = linspace(0, t_end1, data_length1);
tspan2 = linspace(0, t_end2, data_length2);

% Assign coupling strengths
G1 = SyncCouplingAssign(G1, a);
G2 = SyncCouplingAssign(G2, a);

%% Initial conditions
x_mean = 0; 
x_std  = 2;
P      = diag([1, 0, 0]);                    % Projection matrix
X0     = x_mean + x_std*rand(1, numStates*N);

%% Simulate coupled Lorenz systems
% First interval (graph G1)
[X1, t1] = SimulateCoupledSystems(@LorenzOscillator, tspan1, X0, G1, P);

% Use final state from first interval as initial state for second
X0 = X1(end,:);

% Second interval (graph G2)
[X2, t2] = SimulateCoupledSystems(@LorenzOscillator, tspan2, X0, G2, P);

% Concatenate results from both intervals
X = [X1(1:end-1,:); X2];
t = [t1(1:end-1,:); t2 + t1(end)];

%% Visualization settings
state_indices   = 1:numStates;
state_index_all = 1:numStates*N;

colors    = lines(N);
linestyle = {'-','--','-.',':','-','--','-.',':','-','--'};
lw        = [2*ones(1,4) 1.5*ones(1,4) 1 1];
markers   = {'none','none','none','none','*','*','o','o','.','.'};

%% Data storage setup
E         = zeros(N, length(t));   % Synchronization error matrix
cols      = 'ABCDEFGHIJK';         % Excel column labels
filename  = 'sync_data.xlsx';
range_end = length(t) + 1;

%% Plot synchronization errors
figure; hold on; grid on;

for i = 1:N
    % Extract state indices for oscillator i
    slice_i   = (i-1)*numStates + state_indices;
    slice_rem = setdiff(state_index_all, slice_i);

    % Synchronization error (L2 distance from others)
    e = vecnorm(repmat(X(:,slice_i),1,N-1) - X(:,slice_rem), 2, 2).';
    E(i,:) = e;

    % Plot error trajectory
    plot(t, e, ...
        'Color', colors(i,:), ...
        'LineWidth', lw(i), ...
        'LineStyle', linestyle{i}, ...
        'Marker', markers{i}, ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', sprintf('System %d', i));

    % Save to Excel
    range_str = strcat(cols(i+1),'2:',cols(i+1),string(range_end));
    writematrix(e.', filename,'Sheet',1,'Range',range_str);
    writematrix("e"+i, filename,'Sheet',1,'Range',cols(i+1)+"1");
end

xlabel('Time');
ylabel('Synchronization error (L2 norm)');
title('Synchronization of 10 Lorenz Oscillators');
legend show;
hold off;

%% Save time vector to Excel
writematrix(t, filename, 'Sheet', 1, 'Range', strcat(cols(1),'2:',cols(1),string(range_end)));
writematrix('t', filename, 'Sheet', 1, 'Range', 'A1');
