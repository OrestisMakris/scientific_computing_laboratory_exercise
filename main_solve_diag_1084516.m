% Παράδειγμα χρήσης:
clear ;
close all;
A = [2, -1, 0, 0; -1, 2, -1, 0; 0, -1, 2, -1; 0, 0, -1, 2];
B = [1; 2; 3; 4];
D = diag(diag(A));   % Διαγώνιος πίνακας D
L = -tril(A, -1);     % Κάτω τριγωνικό μέρος L
U = -triu(A, 1);      % Άνω τριγωνικό μέρος U
[solution, error_norm, numerical_cost] = solve_tridiagonal_system_1084516(A, B ,D ,  L ,U);

disp('Ιδανίκη Λύση του συστήματος:');
disp(solution);

disp('Νόρμα σφάλματος:');
disp(error_norm);

disp('Αριθμητικό κόστος:');
disp(numerical_cost);

disp('-----------------------------------------------')
disp('-----------------------------------------------')
disp('-----------Διαγώνια κυρίαρχα μητρώα-----------')
disp('-----------------------------------------------')
disp('-----------------------------------------------')
% Πείραμα με τυχαίους διαγώνια κυρίαρχους πίνακες
n_values = [10];
for n = n_values
    A_diag_dom = diag(2866*ones(1, n)) + diag(-ones(1, n-1), 1) + diag(-ones(1, n-1), -1);
    B_diag_dom = randn(n, 1);

    % Διαίρεσε τον πίνακα A σε D, L, U
    D_diag_dom = diag(diag(A_diag_dom));
    L_diag_dom = -tril(A_diag_dom, -1);
    U_diag_dom = -triu(A_diag_dom, 1);

    [solution_diag_dom, error_norm_dom , numerical_cost_diag_dom] = solve_tridiagonal_system_1084516(A_diag_dom, B_diag_dom, D_diag_dom, L_diag_dom, U_diag_dom);

    disp(['Για n = ', num2str(n)]);
    disp(['Αριθμητικό κόστος:', num2str(numerical_cost_diag_dom)]);

    % Υπολογισμός της 2-νόρμας του σφάλματος
    error_norm_diag_dom = norm(solution_diag_dom - A_diag_dom \ B_diag_dom);
    disp(['Νόρμα σφάλματος:', num2str(error_norm_diag_dom)]);

    % Υπολογισμός των νορμών των υπερ/υπό-διαγωνίων
    norm_LD_invU = norm(L_diag_dom * inv(D_diag_dom) * U_diag_dom);
    norm_UD_invL = norm(U_diag_dom * inv(D_diag_dom) * L_diag_dom);

    disp(['Νόρμα του πίνακα LD^(-1)U:', num2str(norm_LD_invU)]);
    disp(['Νόρμα του πίνακα UD^(-1)L:', num2str(norm_UD_invL)]);
end
