function [X,t] = SimulateCoupledSystems(systemDynamics,tSpan,X0,G,P)
    
    A = full(adjacency(G,"weighted")).';

    % Compute the Laplacian matrix from the adjacency matrix
    L = diag(sum(A,2))-A; % In-degree Laplacian
   

    % Number of nodes (oscillators)
    N = G.numnodes;
    % disp(L*ones(N,1))
    stateDim = length(X0)/N; % Number of states per node

    % % Ensure initial condition is properly reshaped
    % X0 = reshape(X0, [], 1); % Convert to column vector

    % Solve the coupled system using ode45
    [t, X] = ode45(@(t,X) coupledDynamics(t,X,systemDynamics,L,P), tSpan, X0);

       
    % Compute the norm of each node's state vector
    % stateNorms = vecnorm(reshape(X, stateDim, N, []), 2, 1);
    % stateNorms = squeeze(stateNorms)'; % Ensure proper dimensions

    % Plot individual state trajectories
    % colors = {'-r';':g';'-.b';'--c';':m';'-y';':r';'-.g';'--b';'c'};
    % 
    % for i = 1:stateDim
    %     figure
    %     grid on
    %     hold on
    %     for j = 1:N
    %         if j>5
    %             p = 1;
    %         else
    %             p = 2;
    %         end
    %         plot(t, X(:,stateDim*(j-1)+i).',colors{j}(:), 'LineWidth',p, 'DisplayName', sprintf('Node %d', j))
    %         legend show
    %     end
    %     xlabel('Time (t)');
    %     ylabel(sprintf('x_%d(t)',i));
    %     title('State Trajectories of Diffusively Coupled Systems');
    %     hold off
    % end

    % Labels and title
    
    
end