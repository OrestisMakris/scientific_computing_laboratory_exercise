function z = preconditioner_my1084516(r, A)
 
    % Εκτέλεση ατελούς παραγοντοποίησης Cholesky
    L = ilu(A);

    % Επίλυση του γραμμικού συστήματος χρησιμοποιώντας ατελή παράγοντα Cholesky
    z = L' \ (L \ r);
end
