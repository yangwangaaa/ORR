function [Hbar, gbar, Abar, bbar, Mbar, mbar] = LTV_MPC_to_QP(x0, xguess, WN, JN, W, J, A, B, r, C, D, h, HN, hN)
%LTV-MPC Converts linear MPC problem with time varying matrices to QP
%   
% This function converts MPC problem in the following form
%
% min_{Dx_i, Du_i} 1/2*Dx_N'W_N*Dx_N + J_N'*Dx_N  sum_{i=0}^{N-1} 1/2*[Dx_i; Du_i]'*W_i*[Dx_i; Du_i] + J_i'*[Dx_i; Du_i]
%             s.t. Dx_0     = x0 - xguess_0,
%                  Dx_{i+1} = A_i*Dx_i + B_i*Du_i + r_i,    i=0,...,N-1,
%                  C_i*Dx_i + D_i*Du_i + h_i <= 0,          i=0,...,N-1,
%                  H_N*Dx_N + h_N <= 0,
%
% where Dx_i=[x_i-xref_i] and Du_i=[u_i-uref_i] to the following QP:
%
% min_{z} 1/2z'*Hbar*z + gbar'*z
%    s.t. Abar*z  = bbar,
%         Mbar*z <= mbar,
% where z = [Dx_0, Du_0, ..., Dx_{N-1}, Du_{N-1}, Dx_N].

% control horizon
N = numel(A);

% Number of states
n = size(A{1}, 1);

% Number of control inputs
m = size(B{1}, 2);

% The dynamics constraints
Abar = kron(eye(N+1), [eye(n) zeros(n,m)]);
Abar = Abar(1:end, 1:end-m);

bbar = zeros((N+1)*n,1);

% The Hessian matrix
Hbar = zeros( (N+1)*n + m*N );
gbar = zeros( (N+1)*n + m*N, 1 );

% Fill in the part of Hbar and gbar corresponding to WN and JN
Hbar((end-n+1):end,(end-n+1):end) = WN;
gbar((end-n+1):end) = JN;

% Fill in the part of bbar corresponding to x0
bbar(1:n) = x0 - xguess(:, 1);
% and r vector
bbar((n+1):end) = r(:);

for i=1:N
    % A_i matrix
    Abar( i*n + (1:n), (i-1)*(n+m) + (1:n))                = -A{i};
    % B_i matrix
    Abar( i*n + (1:n), (i-1)*(n+m) + n + (1:m))            = -B{i};
    
    % H_i
    Hbar( (i-1)*(n+m) + (1:(n+m)), (i-1)*(n+m) + (1:(n+m)) ) = W{i};
    % J_i
    gbar( (i-1)*(n+m) + (1:(n+m)), 1 )                       = J{i};
end

% The inequality constraints
% Compute the the total number of the constraints
n_constr = 0;
for i=1:N
    n_constr = n_constr + size(C{i}, 1);
end
n_constr = n_constr + size(hN, 1);

Mbar = zeros(n_constr, (N+1)*n + N*m);
mbar = zeros(n_constr, 1);

k = 0;
for i=1:N
    % Number of constraints for at ith time step
    nk_constr = size(C{i}, 1);
    if nk_constr < 1
        continue;
    end
    
    Mbar(k + (1:nk_constr), (i-1)*(n+m) +(1:n)) = C{i};
    Mbar(k + (1:nk_constr), (i-1)*(n+m) + n + (1:m)) = D{i}; 
    mbar(k + (1:nk_constr)) = -h{i};
    
    k = k + nk_constr;
end
% Add the final-time path constraint
if size(hN, 1) > 0
    % Number of constraints for at ith time step
    nk_constr = size(hN, 1);
    
    Mbar(k + (1:nk_constr), N*(n+m) +(1:n)) = HN; 
    mbar(k + (1:nk_constr)) = -hN;
end

end

