%% Script έλεγχου και σύγκρισης επιλυτών για (A+u*v')^k * x = b
% Αυτό το script συγκρίνει πέντε διαφορετικές μεθόδους επίλυσης του 
% συστήματος (A+u*v')^k * x = b για διάφορες τιμές του n.
%
% Τελική αξιολόγηση:
% Σε κανονικές συνθήκες τα σφάλματα πρέπει να είναι τάξης O(κ(A+u*v')*eps) 
% (ή μικρότερα), όπου κ(A+u*v') είναι ο δείκτης κατάστασης του συστήματος και 
% eps η μηχανή ακρίβειας.
%
% Το script υπολογίζει επίσης την τιμή κ (χρησιμοποιώντας condest για μεγάλους πίνακες)
% και την θεωρητική τιμή του σφάλματος: theoreticalError = κ(A+u*v') * eps.
% Αυτή η τιμή εκτυπώνεται μαζί με τα πραγματικά σφάλματα για κάθε μέθοδο.

clear; clc;

ks = 5;             % Επιλέγουμε k = 5 (π.χ.)
ns = [500, 1000, 2000];
methodNames = {'SLV1', 'SLV2', 'SLV3', 'SLVpcg', 'SLVnonsym'};
errors = zeros(length(ns), length(methodNames));
times  = zeros(length(ns), length(methodNames));

for idx = 1:length(ns)
    n = ns(idx);
    fprintf('n = %d\n', n);
    
    % Κατασκευή του τριδιαγώνιου, αραιού και ΣΘΟ πίνακα A:
    e = ones(n,1);
    A = spdiags([-e, 4*e, -e], -1:1, n, n);
    
    % Σταθερή αρχικοποίηση
    rng(0);
    u = randn(n,1);
    v = randn(n,1);
    x_true = ones(n,1);

    
    % Κατασκευή του δεξιού μέλους b:
    % b = (A + u*v')^k * x_true, χωρίς σχηματισμό του (A+u*v') ρητά.
    b = x_true;
    for j = 1:ks
        b = A * b + u * (v' * b);
    end
    
    % Υπολογισμός του δείκτη κατάστασης του (A+u*v') χρησιμοποιώντας condest.
    % Για μεγάλους πίνακες είναι προτιμότερο το condest έναντι του cond.
    matCond = condest(A + u*v');
    theoreticalError = matCond * eps;
    fprintf('Condition number of A+u*v'': %e\n', matCond);
    fprintf('Theoretical error bound (κ*eps): %e\n', theoreticalError);
    

    %% Μέθοδος 1: LowRankPowerSLV1_id
    tic;
    x1 = LowRankPowerSLV1_1084516(A, u, v, ks, b);
    times(idx,1) = toc;
    errors(idx,1) = norm(x1 - x_true, inf);
    
    %% Μέθοδος 2: LowRankPowerSLV2_id
    tic;
    x2 = LowRankPowerSLV2_1084516(A, u, v, ks, b);
    times(idx,2) = toc;
    errors(idx,2) = norm(x2 - x_true, inf);
    

    %% Μέθοδος 3: LowRankPowerSLV3_id
    tic;
    x3 = LowRankPowerSLV3_1084516(A, u, v, ks, b);
    times(idx,3) = toc;
    errors(idx,3) = norm(x3 - x_true, inf);
    

    %% Μέθοδος 4: LowRankPowerSLVpcg_id (με χρήση pcg) 
    tic;
    x4 = LowRankPowerSLVpcg_1084516(A, u, v, ks, b);
    times(idx,4) = toc;  

    errors(idx,4) = norm(x4 - x_true, inf); 
    
    %% Μέθοδος 5: LowRankPowerSLVnonsym_id (με GMRES)
    tic;
    x5 = LowRankPowerSLVnonsym_1084516(A, u, v, ks, b);
    times(idx,5) = toc;
    errors(idx,5) = norm(x5 - x_true, inf);
    

    % Εμφάνιση αποτελεσμάτων για το τρέχον n: 
    fprintf('Σφάλματα:\n');
    for m = 1:length(methodNames)
        fprintf('  %s:  %e\n', methodNames{m},  errors(idx, m));
    end 

    fprintf('Χρόνοι (sec):\n');

    for m = 1:length(methodNames)
        fprintf('  %s: %f\n', methodNames{m}, times(idx, m));
    end

    fprintf('---------------------------- ----\n\n');
    
end