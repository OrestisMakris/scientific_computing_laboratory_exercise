% Function using BiCGSTAB for the matrix-free Krylov method
function x = LowRankPowerSLVnonsym_1084516(A, u, v, k, b)
% LowRankPowerSLVnonsym_bicgstab_id
% Solves the system G*x = b, where G = (A+u*v')^k, using a matrix-free
% implementation of the BiCGSTAB method.
%
% This function avoids forming G explicitly by computing its action on vectors
% through repeated multiplication using power_mult.
%
% Inputs:
%   A  : n×n sparse matrix
%   u,v: n×1 vectors
%   k  : positive integer indicating the power
%   b  : right-hand side vector
%
% Output:
%   x  : approximate solution of (A+u*v')^k x = b
%
% The method uses a stopping tolerance with relative residual norm below 10^-6.

    tol = 1e-6;
    maxit =1000;  % maximum iterations (can be adjusted)
    
    % Define the matrix-free operator as a function handle.
    Ghandle = @(x) power_mult(x, A, u, v, k);
    
    % Use BiCGSTAB to solve G*x = b.
    [x, flag, relres] = bicgstab(Ghandle, b, tol, maxit);
    if flag ~= 0
        warning('BiCGSTAB did not converge: flag = %d, relres = %e', flag, relres);
    end
end

%% Helper function for matrix-free multiplication with G = (A+u*v')^k.
function y = power_mult(x, A, u, v, k)
    % Computes y = (A+u*v')^k * x without explicitly forming the matrix (A+u*v')^k.
    y = x;
    for i = 1:k
        y = A * y + u * (v' * y);
    end
end