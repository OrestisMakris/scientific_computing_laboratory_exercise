clear; clc;

%% Test for small n = 4
n = 4;
% Construct A as a Toeplitz matrix with first column [4; -1; zeros(n-2,1)]
c = [4; -1; zeros(n-2, 1)];
r = [4, -1, zeros(1, n-2)];
A = toeplitz(c, r);

rng(0); % Set the seed for reproducibility
u = randn(n,1);
v = randn(n,1);
x = ones(n,1);
k = 3;   % Note: 1 < k < n

% Direct computation of (A+u*v')^k * x:
A_direct = (A + u*v')^k;
b_direct = A_direct * x;

% Compute U,V using the low-rank method:
[U, V] = LowRankPower_1084516(A, u, v, k);

% Reconstruct (A+u*v')^k using the low-rank identity:
A_reconstructed = A^k + U*V';
b_lowrank = A_reconstructed * x;

% Compute the infinity norm errors:
err_b = norm(b_direct - b_lowrank, inf);
err_A = norm(A_direct - A_reconstructed, inf);

fprintf('For n = %d:\n', n);
fprintf('Infinity norm of error in b: %e\n', err_b);
fprintf('Infinity norm of error in A: %e\n\n', err_A);

%% Timing experiments for larger n
ns = [500, 1000, 2000];
timings = struct('n', [], 'direct', [], 'lowrank', [], 'errA', [], 'errb', []);

for idx = 1:length(ns)
    n = ns(idx);
    % Create Toeplitz matrix A (n-by-n)
    c = [4; -1; zeros(n-2, 1)];
    r = [4, -1, zeros(1, n-2)];
    A = toeplitz(c, r);
    
    rng(0); % Ensure same random data for each run.
    u = randn(n,1);
    v = randn(n,1);
    x = ones(n,1);
    k = 3; % keeping k fixed
    
    % Time direct computation
    tic;
    A_direct = (A + u*v')^k;
    b_direct = A_direct * x;
    time_direct = toc;
    
    % Time low-rank computation
    tic;
    [U, V] = LowRankPower_1084516(A, u, v, k);
    A_reconstructed = A^k + U*V';
    b_lowrank = A_reconstructed * x;
    time_lowrank = toc;
    
    % Compute errors
    err_A = norm(A_direct - A_reconstructed, inf);
    err_b = norm(b_direct - b_lowrank, inf);
    
    % Save results
    timings(idx).n = n;
    timings(idx).direct = time_direct;
    timings(idx).lowrank = time_lowrank;
    timings(idx).errA = err_A;
    timings(idx).errb = err_b;
    
    fprintf('For n = %d:\n', n);
    fprintf('  Direct method:    %f seconds\n', time_direct);
    fprintf('  Low-rank method:  %f seconds\n', time_lowrank);
    fprintf('  Infinity norm error in A: %e\n', err_A);
    fprintf('  Infinity norm error in b: %e\n\n', err_b);
end