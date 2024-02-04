% Example usage:
A = [2, -1, 0, 0; -1, 2, -1, 0; 0, -1, 2, -1; 0, 0, -1, 2]
B = [1; 2; 3; 4];

solution = solve_tridiagonal_system(A, B);
disp('Solution of the system:');
disp(solution);

x_direct= A \ B 
error_norm = norm(solution - x_direct)

