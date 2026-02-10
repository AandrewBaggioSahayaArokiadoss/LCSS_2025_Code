% Synchronization of a dynamical network of Lorenz oscillators
% ------------------------------------------------------------
% Purpose:
%   Simulate a network of N identical Lorenz oscillators coupled via a directed graph.
%   Coupling strengths are assigned via SyncCouplingAssign to guarantee synchronization.
%   After simulation, compute pairwise synchronization error (state-space L2 distance)
%   between each oscillator and its "next" in a ring (1–2, 2–3, …, N–1–N, N–1).
%   Plot the error vs time with distinct line/marker styles for each pair.
%
% Note:
%   P (projection matrix) determines which state variables (e.g. x-component) are coupled.
%   Here P = diag([1 0 0]) → only the first coordinate is used in coupling.

clc; clear; close all;

%% === System parameters (Lorenz oscillator) ===
sigma = 10;
rho   = 25;
beta  = 8/3;

% Analytical coupling strength derived (sufficient for synchronization)
a = -sigma + (beta*(beta+1)*(rho+sigma)^2) / (16*(beta-1));

%% === Define connectivity digraph ===
N = 10;       % Number of oscillators
density = 0.1;
G = generateOneRootSCC(N, density);
num_states = 3;        % Dimension of each oscillator's state

%% === Assign coupling strengths ===
G = SyncCouplingAssign(G, 1.1*a);

figure;
plot(G, 'EdgeLabel', G.Edges.Weight, 'Layout', 'circle');
title('Coupling digraph with assigned weights');

%% === Projection matrix for coupling ===
P = diag([1, 0, 0]);   % couple only first coordinate (e.g. x-variable)

%% === Simulation settings ===

data_length = 10;
t_end = 50;
tspan = linspace(0, t_end, data_length);

%% === Initial conditions ===
x_mean = 0;
x_std  = 10^-6;
X0 = x_mean + x_std * rand(num_states * N,1);

%% === Simulate coupled Lorenz oscillators ===
[X, t] = SimulateCoupledSystems(@LorenzOscillator,tspan, X0, G, P,a);

% INPUT:
% t             → d×1 vector of time points
% X             → d×(N*num_states) state data
% N             → number of oscillators
% num_states    → number of states per oscillator

X = X.';

X_pdist = zeros(N,length(t));

A = adjacency(digraph(1:N-1,2:N));
L_pair = diag(sum(A,2)) - A;
L_kron = kron(L_pair,eye(num_states));

% compute pairwise distances and sum for each time step
for ti = 1:length(t)
    X_pdist(:,ti) = sqrt(sum(reshape(L_kron*X(:,ti),num_states,N).^2,1)).';
end

% plot
figure;
plot(t,X_pdist, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Sum of All Pairwise Distances');
title('Total Pairwise Distance vs Time');
grid on;
