function x = solve_tridiagonal_system(A, B)
    % Solve tridiagonal system Ax = 
    % Split the matrix A into D, L, U
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, 1);
    UU = D - L - U  ;
    inv_D = inv(D);
    YUU = U*inv_D*L

    M = (D + L + U) * inv_D * (D - L - U);

    b_prime = (D + L + U) * inv_D * B;

% Λύνουμε το σύστημα a * x = b'
x = M \ b_prime;

end


