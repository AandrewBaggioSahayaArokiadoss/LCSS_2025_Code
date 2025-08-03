function dxdt = LorenzOscillator(t,X)
    sigma = 10;
    rho = 28;
    beta = 8/3;
    dxdt = [sigma * (X(2) - X(1)); X(1) * (rho - X(3)) - X(2); X(1) * X(2) - beta * X(3)];
end
