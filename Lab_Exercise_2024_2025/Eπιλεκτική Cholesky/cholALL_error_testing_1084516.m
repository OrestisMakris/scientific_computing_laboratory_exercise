n = 1000;
rng('default');
d = rand(n,1);
u = randn(n,1);
A = diag(d) + u*u';

% Υπολογισμός των διαγώνιων και υποδιαγώνιων τιμών με τη συνάρτηση cholALL_ID
[ld, ls] = cholALL_1084516(u,d);

% Υπολογισμός της παραγοντοποίησης Cholesky μέσω της MATLAB
L_builtin = chol(A, 'lower');

% Υπολογισμός της μέγιστης απόκλισης για τις διαγώνιες και υποδιαγώνιες τιμές
max_diff_diag = max(abs(ld - diag(L_builtin)));
max_diff_subdiag = max(abs(ls(2:end) - diag(L_builtin, -1)));

% Εκτύπωση αποτελεσμάτων
fprintf('Maximal difference in diagonal: %e\n', max_diff_diag);
fprintf('Maximal difference in subdiagonal: %e\n', max_diff_subdiag);
if max_diff_diag < n*eps && max_diff_subdiag < n*eps
    fprintf('Ο υπολογισμός είναι ακριβής (στο επίπεδο της μηχανικής ακρίβειας).\n');
else
    fprintf('Υπάρχει σφάλμα στον υπολογισμό των στοιχείων του L.\n');
end

% Υπολογισμός όλων των διαγώνιων και υποδιαγώνιων τιμών
max_diff_diag_full = max(abs(diag(L_builtin) - ld));
max_diff_subdiag_full = max(abs(diag(L_builtin, -1) - ls(2:end)));

% Εκτύπωση αποτελεσμάτων για όλες τις διαγώνιες και υποδιαγώνιες τιμές
fprintf('Maximal difference in all diagonal elements: %e\n', max_diff_diag_full);
fprintf('Maximal difference in all subdiagonal elements: %e\n', max_diff_subdiag_full);

% Εξαγωγή της μέγιστης απόκλισης για κάθε τιμή (από j=1 έως j=n για την διαγώνιο)
max_diff_diag_all = max(abs(diag(L_builtin) - ld));
% Εξαγωγή της μέγιστης απόκλισης για την υποδιαγώνιο (από j=2 έως j=n-1)
max_diff_subdiag_all = max(abs(diag(L_builtin, -1) - ls(2:end)));

% Εκτύπωση της απόκλισης για όλες τις διαγώνιες και υποδιαγώνιες τιμές
fprintf('Maximal difference in diagonal elements (all): %e\n', max_diff_diag_all);
fprintf('Maximal difference in subdiagonal elements (all): %e\n', max_diff_subdiag_all);

% Έλεγχος αν το σφάλμα είναι μικρότερο από το O(neps)
if max_diff_diag_all < n*eps && max_diff_subdiag_all < n*eps
    fprintf('Ο υπολογισμός των στοιχείων του L είναι σωστός (στο επίπεδο ακρίβειας O(neps)).\n');
else
    fprintf('Υπάρχει σφάλμα στον υπολογισμό των στοιχείων του L.\n');
end
