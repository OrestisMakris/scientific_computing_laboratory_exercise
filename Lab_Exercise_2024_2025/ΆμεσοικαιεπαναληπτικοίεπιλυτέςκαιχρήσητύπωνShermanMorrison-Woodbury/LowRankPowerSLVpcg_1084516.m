%% Function using PCG (with incomplete Cholesky preconditioner) for solves with A.
function x = LowRankPowerSLVpcg_1084516(A, u, v, k, b)
% LowRankPowerSLVpcg_id
% Solves (A+u*v')^k * x = b using an iterative approach where at each step,
% the linear system with A is solved using the preconditioned conjugate gradient method (pcg)
% with an incomplete Cholesky preconditioner. The tolerance is set so that the relative
% residual norm stays below 10^-6.
%
% Inputs:
%   A  : sparse symmetric positive definite tridiagonal matrix (n×n)
%   u,v: column vectors (n×1)
%   k  : positive integer (power)
%   b  : right-hand side vector (n×1)
%
% Output:
%   x  : approximate solution of (A+u*v')^k * x = b

    tol = 1e-6;
    maxit = min(size(A,1), 100);

    % Compute an incomplete Cholesky factorization of A for preconditioning.
    L = ichol(A);
    
    % Pre-solve A * x = u efficiently with preconditioned CG.
    [A_inv_u, flag] = pcg(A, u, tol, maxit, L, L');
    if flag ~= 0
        warning('pcg did not converge for solving A\\u.');
    end
    denom = 1 + v' * A_inv_u;
    
    y = b;
    for i = 1:k
        [z, flag] = pcg(A, y, tol, maxit, L, L');
        if flag ~= 0
            warning('pcg did not converge for solving A\\y on iteration %d.', i);
        end
        correction = (v' * z) / denom;
        y = z - A_inv_u * correction;
    end
    x = y;
end