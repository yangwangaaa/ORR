function x_new = rk4_int(f, x, u, Ts)
    % Runge-Kutta 4 integration
    k1 = f(x,         u);
    k2 = f(x+Ts/2*k1, u);
    k3 = f(x+Ts/2*k2, u);
    k4 = f(x+Ts*k3,   u);
    x_new = x + Ts/6*(k1+2*k2+2*k3+k4);
end