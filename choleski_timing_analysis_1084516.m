function choleski_timing_analysis_1084516
    clear;
    close all;
    
    % Ορισμός των τιμών του n και n1
    n_values = 100:100:2000;
    n1_values = 150:100:1550;
    
    % Πλήθος των πινάκων
    num_matrices = numel(n_values);
    
    % Δημιουργία κελιού για την αποθήκευση των πινάκων
    matrices_cell_A = cell(1, num_matrices);
    
    % Δημιουργία κελιού για τους χρόνους Cholesky
    choleski_times = zeros(size(n_values));
    
    % Λούπα για τη δημιουργία των τυχαίων συμμετρικών και θετικά ορισμένων πινάκων
    for i = 1:num_matrices
        n = n_values(i);
        
        % Δημιουργία τυχαίου πίνακα μεγέθους n x n
        A = randn(n, n);
        
        % Κάντε τον πίνακα συμμετρικό
        A = 0.5 * (A + A');
        
        % Κάντε τον πίνακα θετικά ορισμένο
        A = A + n * eye(n);
        
        % Αποθήκευση του πίνακα στο κελί
        matrices_cell_A{i} = A;
    end
    
    % Λούπα για τον υπολογισμό των χρόνων Cholesky
    for i = 1:num_matrices
        choleski_times(i) = timeit(@() chol(matrices_cell_A{i}));
    end
    
    p = polyfit(n_values,  choleski_times, 3);
    Tchol = @(n) polyval(p, n);
    
    
    % Πρόβλεψη χρόνων για τις τιμές n_values και n1_values
    times_predicted = Tchol(n_values);
    times_predicted1 = Tchol(n1_values);
    
    % Υπολογισμός απόλυτης απόκλισης (ακρίβειας) για n = [100:100:2000]
    actual_times_n = choleski_times;
    error_n = abs(actual_times_n - times_predicted);
    
    % Υπολογισμός απόλυτης απόκλισης (ακρίβειας) για n = [150:100:1550]
    actual_times_n1 = zeros(size(n1_values));
    for i = 1:numel(n1_values)
        actual_times_n1(i) = timeit(@() chol(randn(n1_values(i), n1_values(i)) + n1_values(i) * eye(n1_values(i))));
    end
    error_n1 = abs(actual_times_n1 - times_predicted1);
    
    % Εμφάνιση της ακρίβειας για κάθε πρόβλεψη
    fprintf('Ακρίβεια για n = [100:100:2000]: %.4f\n', mean(error_n));
    fprintf('Ακρίβεια για n = [150:100:1550]: %.4f\n', mean(error_n1));
    
    % Οπτικοποίηση αποτελεσμάτων
    figure;
    hold on;
    
    % Χρονομετρημένοι χρόνοι για n = [100:100:2000]
    plot(n_values, actual_times_n, 'o', 'DisplayName', 'Μετρημένοι Χρόνοι (n = [100:100:2000])');
    
    % Χρονομετρημένοι χρόνοι για n = [150:100:1550]
    plot(n1_values, actual_times_n1, 'o', 'DisplayName', 'Μετρημένοι Χρόνοι (n = [150:100:1550])');
    
    % Προβλεπόμενοι χρόνοι από τη συνάρτηση πρόβλεψης
    plot(n_values, times_predicted, 'DisplayName', 'Προβλεπόμενοι Χρόνοι (n = [100:100:2000])', 'LineWidth', 2);
    
    % Προβλεπόμενοι χρόνοι για n = [150:100:1550]
    plot(n1_values, times_predicted1, 'DisplayName', 'Προβλεπόμενοι Χρόνοι (n = [150:100:1550])', 'LineWidth', 2);
    
    xlabel('Μέγεθος Πίνακα (n)');
    ylabel('Cholesky Χρόνος');
    title('Ακρίβεια Πρόβλεψης Χρόνου Cholesky');
    legend('show');
    grid on;
    hold off;
    
    % Πολυώνυμο 2ου βαθμού
    p2 = polyfit(n_values,  choleski_times , 2);
    Tchol2 = @(n) polyval(p2, n);
    times_predicted2 = Tchol2(n_values);
    times_predicted2_1 = Tchol2(n1_values);
    
    % Πολυώνυμο 4ου βαθμού
    p4 = polyfit(n_values,  choleski_times , 4);
    Tchol4 = @(n) polyval(p4, n);
    times_predicted4 = Tchol4(n_values);
    times_predicted4_1 = Tchol4(n1_values);
    
    % Υπολογισμός απόλυτης απόκλισης (ακρίβειας) για n = [100:100:2000] για πολυώνυμα 2ου και 4ου βαθμού
    error_n2 = abs(actual_times_n - times_predicted2);
    error_n4 = abs(actual_times_n - times_predicted4);
    
    % Υπολογισμός απόλυτης απόκλισης (ακρίβειας) για n = [150:100:1550] για πολυώνυμα 2ου και 4ου βαθμού
    error_n2_1 = abs(actual_times_n1 - times_predicted2_1);
    error_n4_1 = abs(actual_times_n1 - times_predicted4_1);
    
    % Εμφάνιση της ακρίβειας για κάθε πρόβλεψη
    fprintf('Ακρίβεια για n = [100:100:2000] (Τετραγωνικό): %.4f\n', mean(error_n2));
    fprintf('Ακρίβεια για n = [100:100:2000] (Τεταρτοβαθμιαίο): %.4f\n', mean(error_n4));
    fprintf('Ακρίβεια για n = [150:100:1550] (Τετραγωνικό): %.4f\n', mean(error_n2_1));
    fprintf('Ακρίβεια για n = [150:100:1550] (Τεταρτοβαθμιαίο): %.4f\n', mean(error_n4_1));
    
    % Οπτικοποίηση αποτελεσμάτων
    figure;
    
    % Πρώτο subplot για το πολυώνυμο 2ου βαθμού
    subplot(2, 2, 1);
    hold on;
    plot(n_values, actual_times_n, 'o', 'DisplayName', 'Μετρημένοι Χρόνοι (n = [100:100:2000])');
    plot(n_values, times_predicted2, 'DisplayName', 'Προβλεπόμενοι Χρόνοι (Τετραγωνικό)', 'LineWidth', 2);
    xlabel('Μέγεθος Πίνακα (n)');
    ylabel('Cholesky Χρόνος');
    title('Ακρίβεια Πρόβλεψης (Τετραγωνικό) - n = [100:100:2000]');
    legend('show');
    grid on;
    hold off;
    
    % Δεύτερο subplot για το πολυώνυμο 2ου βαθμού (n = [150:100:1550])
    subplot(2, 2, 2);
    hold on;
    plot(n1_values, actual_times_n1, 'o', 'DisplayName', 'Μετρημένοι Χρόνοι (n = [150:100:1550])');
    plot(n1_values, times_predicted2_1, 'DisplayName', 'Προβλεπόμενοι Χρόνοι (Τετραγωνικό)', 'LineWidth', 2);
    xlabel('Μέγεθος Πίνακα (n)');
    ylabel('Cholesky Χρόνος');
    title('Ακρίβεια Πρόβλεψης (Τετραγωνικό) - n = [150:100:1550]');
    legend('show');
    grid on;
    hold off;
    
    % Τρίτο subplot για το πολυώνυμο 4ου βαθμού
    subplot(2, 2, 3);
    hold on;
    plot(n_values, actual_times_n, 'o', 'DisplayName', 'Μετρημένοι Χρόνοι (n = [100:100:2000])');
    plot(n_values, times_predicted4, 'DisplayName', 'Προβλεπόμενοι Χρόνοι (Τεταρτοβαθμιαίο)', 'LineWidth', 2);
    xlabel('Μέγεθος Πίνακα (n)');
    ylabel('Cholesky Χρόνος');
    title('Ακρίβεια Πρόβλεψης (Τεταρτοβαθμιαίο) - n = [100:100:2000]');
    legend('show');
    grid on;
    hold off;
    
    % Τέταρτο subplot για το πολυώνυμο 4ου βαθμού (n = [150:100:1550])
    subplot(2, 2, 4);
    hold on;
    plot(n1_values, actual_times_n1, 'o', 'DisplayName', 'Μετρημένοι Χρόνοι (n = [150:100:1550])');
    plot(n1_values, times_predicted4_1, 'DisplayName', 'Προβλεπόμενοι Χρόνοι (Τεταρτοβαθμιαίο)', 'LineWidth', 2);
    xlabel('Μέγεθος Πίνακα (n)');
    ylabel('Cholesky Χρόνος');
    title('Ακρίβεια Πρόβλεψης (Τεταρτοβαθμιαίο) - n = [150:100:1550]');
    legend('show');
    grid on;
    hold off;
end 