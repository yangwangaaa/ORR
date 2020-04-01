function x_new = euler_int(f, x, u, Ts)
    % Explicit Euler
    x_new = x + Ts*f(x, u);
end