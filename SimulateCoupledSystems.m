% The system takes in the ODE of an oscillator, the connectivity digraph 
% and the P matrix (in paper 10.1109/TCSI.2015.2395632)
% It returns the state vector for all the dynamical systems in the network

function [X,t] = SimulateCoupledSystems(systemDynamics,tSpan,X0,G,P)   

    % Number of nodes (oscillators)
    N = G.numnodes;

    A = full(adjacency(G,"weighted")).';

    % Compute the Laplacian matrix from the adjacency matrix
    L = diag(sum(A,2))-A; % In-degree Laplacian

    % Solve the coupled system using ode45
    [t, X] = ode45(@(t,X) coupledDynamics(t,X,systemDynamics,L,P), tSpan, X0);
end
