function result = ttt_myid_1084516(A, B, varargin)
    % Έλεγχος για τον αριθμό των ορισμάτων
    if nargin < 2
        error('Λάθος αριθμός ορισμάτων. Χρειάζονται τουλάχιστον δύο τανυστές.');
    end
    
    % Έλεγχος για τον τύπο των εισόδων
    if ~isnumeric(A) || ~isnumeric(B)
        error('Οι τανυστές πρέπει να είναι αριθμητικοί πίνακες.');
    end
    
    % Έλεγχος για τις διαστάσεις των τανυστών
    if ~isequal(size(A), size(B))
        % error('Οι τανυστές πρέπει να έχουν ίδιες διαστάσεις.');
    end
    
    % Έλεγχος για το είδος της πράξης (εσωτερικό ή εξωτερικό γινόμενο)
    if nargin == 3 && strcmpi(varargin{1}, 'all')
        % Εσωτερικό γινόμενο
        result = sum(A(:) .* B(:));
    else
        % Αναδιάταξη και υπολογισμός του εξωτερικού γινομένου
        result = reshape(A, [], 1) * reshape(B, 1, []);
        result = reshape(result, [size(A), size(B)]); % Μέγεθος εξόδου [3, 3, 3, 3, 3, 3]
    end
end