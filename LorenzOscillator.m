% LorenzOscillator defines the Lorenz system ODE.
%
%   dxdt = LorenzOscillator(t,X)
%
%   Inputs:
%     t – time (unused, included for ode45 compatibility)
%     X – state vector [x; y; z]
%
%   Outputs:
%     dxdt – time derivative of state vector

function dxdt = LorenzOscillator(~,X)

    %% Parameters
    sigma = 10;
    rho   = 28;
    beta  = 8/3;

    %% Dynamics
    dxdt = [
        sigma * (X(2) - X(1));           % dx/dt
        X(1) * (rho - X(3)) - X(2);      % dy/dt
        X(1) * X(2) - beta * X(3)        % dz/dt
    ];
end
