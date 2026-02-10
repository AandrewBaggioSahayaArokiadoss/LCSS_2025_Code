% This function returns dX/dt for the entire network by evaluating the
% state of each system

function dXdt = CoupledDynamics(t,X,systemDynamics,L,P,a)

N = length(L);
stateDim = size(P,1);
interactions = sum(kron(L,P)*X,2); % Diffusive coupling term
couplingTerm = reshape(interactions,stateDim,N);
X = reshape(X,stateDim,N);

% disp(couplingTerm-couplingTerm2)
dXdt = zeros(stateDim,N);
for i = 1:N
    dXdt(:, i) = systemDynamics(t, X(:, i))-couplingTerm(:,i)-a* P* X(:, i);
end
dXdt = dXdt(:); % Flatten back to column vector
end
