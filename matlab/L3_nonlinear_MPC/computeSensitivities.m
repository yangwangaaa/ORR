function [WN, JN, W, J, A, B, r, C, D, h, HN, hN] = computeSensitivities(x_guess, u_guess, x_ref, u_ref, Ws, WNs, Js, JNs, fs, As, Bs, hs, Cs, Ds, hNs, HNs)
    N = size(u_ref, 2);
    n = size(x_ref, 1);
    
    A = cell(N, 1);
    B = cell(N, 1);
    W = cell(N, 1);
    J = cell(N, 1);
    C = cell(N, 1);
    D = cell(N, 1);
    h = cell(N, 1);
    
    r = zeros(n, N);
    
    for i=1:N
        % Approximation of the dynamics
        A{i} = As(x_guess(:, i), u_guess(:, i));
        B{i} = Bs(x_guess(:, i), u_guess(:, i));
        r(:, i) = fs(x_guess(:, i), u_guess(:, i)) - x_guess(:, i+1);
        
        % Gauss-Newton Hessian approximation
        W{i} = Ws(x_guess(:, i), u_guess(:, i), x_ref(:, i), u_ref(:, i));
        J{i} = Js(x_guess(:, i), u_guess(:, i), x_ref(:, i), u_ref(:, i));
        
        % Approximation of the path constraints
        C{i} = Cs(x_guess(:, i), u_guess(:, i));
        D{i} = Ds(x_guess(:, i), u_guess(:, i));
        h{i} = hs(x_guess(:, i), u_guess(:, i));
    end

    WN = WNs(x_guess(:, N+1), x_ref(:, N+1));
    JN = JNs(x_guess(:, N+1), x_ref(:, N+1));
    HN = HNs(x_guess(:, N+1));
    hN = hNs(x_guess(:, N+1));
end