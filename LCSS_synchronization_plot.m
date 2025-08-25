% Creates the plot for the synchronization of a dynamical network of
% Lorenz oscillators. Also saves synchronization data into an Excel sheet
% named "sync_data.xlsx".

clc;
clear;
close all;

%% Parameters for the Lorenz oscillator
sigma = 10;        % Sigma parameter
rho   = 25;        % Rho parameter
beta  = 8/3;       % Beta parameter

% Coupling strength "a"
a = -sigma + (beta*(beta+1)*(rho+sigma)^2) / (16*(beta-1));

%% Define the directed network (graph)
tail1 = [1 2 3 3 4 8 8 7 8 6 7 9 10 5 5 7];
head1 = [2 3 1 4 1 1 3 3 6 7 8 10 5 9 8 9];
G1    = digraph(tail1,head1);

N1         = G1.numnodes;   % Number of oscillators (nodes)
numStates1 = 3;            % Dimension of each oscillator state

tail2 = [1 2 2 3 3 4 1 4 2 3 4 5 5 6 7 8 8 8 7 9 10];
head2 = [2 1 3 1 4 1 5 5 6 6 7 6 7 8 5 6 7 9 9 10 9];
G2    = digraph(tail2,head2);

N         = 10;   % Number of oscillators (nodes)
numStates = 3;            % Dimension of each oscillator state

%% Simulation settings
data_length1 = 25;
data_length2 = 25;
t_start1     = 0;
t_end1       = 3;
t_start2     = t_end1;
t_end2       = 6;
tspan1       = linspace(0,t_end1,data_length1);
tspan2       = linspace(0,t_end2,data_length2);

% Assign coupling strengths
G1 = SyncCouplingAssign(G1,a);
G2 = SyncCouplingAssign(G2,a);

%% Mean and std for initial conditions
x_mean = 10; 
x_std  = 5;

P  = diag([1,0,0]);                        % Projection matrix
X0 = x_mean + x_std*rand(1,numStates*N);   % Gaussian random initial conditions

% Simulate coupled Lorenz systems
[X1,t1] = SimulateCoupledSystems(@LorenzOscillator,tspan1,X0,G1,P);

X0 = X1(end,:);

[X2,t2] = SimulateCoupledSystems(@LorenzOscillator,tspan2,X0,G2,P);

% X = [X1(1:end-1,:);X2];
% t = [t1(1:end-1,:);t2+t1(end)];

X = X2;
t = t2;

size(X)
size(t)

%% Visualization settings
state_indices      = 1:numStates;
state_index_full = 1:numStates*N;

colors    = lines(N);
linestyle = {'-','--','-.',':','-','--','-.',':','-','--'};
lw        = [2*ones(1,4) 1.5*ones(1,4) 1 1];
markers   = {'none','none','none','none','*','*','o','o','.','.'};

%% Data storage setup
E          = zeros(N,size(t,1));             % Pairwise error storage
cols       = 'ABCDEFGHIJK';                  % Excel column labels
range_end  = size(t,1)+1;
filename   = 'sync_data.xlsx';

%% Plot synchronization errors
figure;
hold on; grid on;

for i = 1:N
    % Indices for this oscillator
    slice_i   = (i-1)*numStates + state_indices;
    slice_rem = setdiff(state_index_full,slice_i);

    % Compute synchronization error
    e = vecnorm(repmat(X(:,slice_i),1,N-1) - X(:,slice_rem),2,2).';
    E(i,:) = e;

    % Plot trajectory of error norm
    plot(t.', e, ...
        'Color',colors(i,:), ...
        'LineWidth',lw(i), ...
        'LineStyle',linestyle{i}(:), ...
        'Marker',markers{i}(:), ...
        'MarkerFaceColor','none', ...
        'DisplayName', sprintf('System %d',i));

    % Write error data into Excel
    range_str = strcat(cols(i+1),'2:',cols(i+1),string(range_end));
    writematrix(e.', filename,'Sheet',1,'Range',range_str);
    writematrix(strcat('e',string(i)), filename,'Sheet',1,'Range',strcat(cols(i+1),'1'));
end

xlabel('Time');
ylabel('Sum of squared pairwise distances');
title('Synchronization of 10 Lorenz Oscillators');
legend show;

hold off;

%% Save time vector to Excel
range_str = strcat(cols(1),'2:',cols(1),string(range_end));
writematrix(t, filename,'Sheet',1,'Range',range_str);
writematrix('t', filename,'Sheet',1,'Range','A1');
