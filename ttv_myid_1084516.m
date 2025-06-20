function Y = ttv_myid_1084516(X, V, N)
    % Έλεγχος εάν η καθορισμένη διάσταση N είναι έγκυρη
    if N < 1 || N > ndims(X)
        error('Μη έγκυρη διάσταση που καθορίζεται.');
    end
    
    % Λήψη του μεγέθους του X
    sz = size(X);
    
    % Έλεγχος εάν το μέγεθος του V είναι συμβατό με την καθορισμένη διάσταση
    if numel(V) ~= sz(N)
        error('Αναντιστοιχία μεγέθους μεταξύ της διάστασης N και του μεγέθους του V.');
    end
    
    % Αναδιάταξη των διαστάσεων του X για να φέρει την καθορισμένη διάσταση στην αρχή
    perm_order = [N, 1:N-1, N+1:ndims(X)];
    X_permuted = permute(X, perm_order);
    
    % Αναδιαμόρφωση του X σε ένα 2D πίνακα για τον πολλαπλασιασμό
    X_reshaped = reshape(X_permuted, sz(N), []);
    
    % Εκτέλεση του πολλαπλασιασμού (βεβαιωθείτε ότι οι διαστάσεις είναι σωστές)
    result_2d = X_reshaped.' * V;  % Τροποποίηση για να ταιριάζουν οι διαστάσεις
    
    % Επαναδιαμόρφωση του αποτελέσματος σε ένα 2D πίνακα χωρίς διαστάσεις singleton
    Y = reshape(result_2d.', size(X, setdiff(1:ndims(X), N)));
end
