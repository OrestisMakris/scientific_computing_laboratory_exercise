clear;
close all;
n_values = 100:100:2000;
n1_values= 150:100:1550;
num_matrices = numel(n_values);
matrices_cell_A = cell(1, num_matrices);
choleski_times = zeros(size(n_values));

for i = 1:num_matrices
    n = n_values(i);
    
    % Δημιουργία τυχαίου πίνακα μεγέθους n x n
    A = randn(n,n);
    
    % Κάντε τον πίνακα συμμετρικό
    A = 0.5 * (A + A');
    
    % Κάντε τον πίνακα θετικά ορισμένο
    A = A + n * eye(n);
    
    % Αποθήκευση του πίνακα στο κελί
    %isposdef = all(eig(A) > 0)
    matrices_cell_A{i} = A;
end

for i = 1:num_matrices
    choleski_times(i) = timeit(@() chol(matrices_cell_A{i}));
end

p = polyfit(n_values,  choleski_times, 3);
Tchol = @(n) polyval(p, n);

times_predicted = Tchol(n_values);
times_predicted1 = Tchol(n1_values);

