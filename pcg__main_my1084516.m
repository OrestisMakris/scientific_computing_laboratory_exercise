% Ορισμός των παραμέτρων
clear
close all;
load('1138_bus.mat');

% Έλεγχος αν το 'A' είναι πεδίο στη δομή που φορτώθηκε
if isfield(Problem, 'A')
    % Πρόσβαση στον πίνακα 'A' από τη δομή
    A = Problem.A;
    disp('Ο πίνακας A φορτώθηκε με επιτυχία.');
else
    disp('Σφάλμα: Ο πίνακας A δεν βρέθηκε στο αρχείο MAT.');
end

% Ορισμός διανύσματος δεξιάς πλευράς b
b = A * ones(1138, 1);

% Τολεράντα για σύγκλιση
tol = 1e-6;

% Μέγιστος αριθμός επαναλήψεων
maxit = 4 * 1138;

% Ταυτοτικός πίνακας ως απλός προϋπολογιστής M1
M1 = eye(size(A));

% Δεν έχει καθοριστεί δεύτερος προϋπολογιστής (M2)
M2 = [];

% Αρχική υπόθεση για τη λύση
x0 = zeros(size(b));

% Κανονικό PCG
tic;
[x_standard, flag_standard, relres_standard, iter_standard, resvec_standard] = pcg(A, b, tol, maxit);
time_standard = toc;

% Καθορισμός της ακριβούς λύσης (για παρακολούθηση σφάλματος)
xsol = x_standard;

% Τροποποιημένο PCG με προϋπολογιστή
tic;
[x, flag, relres, iter, resvec, errvec] = pcg_my1084516(A, b, tol, maxit, M1, M2, x0, xsol);
time_myid = toc;

% Εμφάνιση των αποτελεσμάτων
disp(['Κανονικό PCG:']);
disp(['  Flag: ', num2str(flag_standard)]);
disp(['  Σχετικό Υπόλοιπο: ', num2str(relres_standard)]);
disp(['  Αριθμός Επαναλήψεων: ', num2str(iter_standard)]);
disp(['  Χρόνος Εκτέλεσης: ', num2str(time_standard), ' δευτερόλεπτα']);
disp(' ');

disp(['Τροποποιημένο PCG (με προϋπολογιστή):']);
disp(['  Flag: ', num2str(flag)]);
disp(['  Σχετικό Υπόλοιπο: ', num2str(relres)]);
disp(['  Αριθμός Επαναλήψεων: ', num2str(iter)]);
disp(['  Χρόνος Εκτέλεσης: ', num2str(time_myid), ' δευτερόλεπτα']);
disp(' ');

% Σχεδιασμός ιστορικού σύγκλισης
figure;
semilogy(0:iter_standard, resvec_standard, '-o', 'DisplayName', 'Κανονικό PCG');
hold on;
semilogy(0:iter, resvec, '-o', 'DisplayName', 'Τροποποιημένο PCG');
xlabel('Επανάληψη');
ylabel('Σχετικό Υπόλοιπο');
title('Ιστορικό Σύγκλισης');
legend('show');
grid on;
hold off;

% Σχεδιασμός ιστορικού σφάλματος
figure;
semilogy(0:iter, errvec, '-o', 'DisplayName', 'Α-νόρμα Σφάλματος');
xlabel('Επανάληψη');
ylabel('Μέγεθος');
title('Ιστορικό Σφάλματος');
legend('show');
grid on;

% Εμφάνιση του τελικού σχετικού υπολοίπου και σφάλματος
disp(['Τελικό Σχετικό Υπόλοιπο (Κανονικό PCG): ', num2str(relres_standard)]);
disp(['Τελικό Σχετικό Υπόλοιπο (Τροποποιημένο PCG): ', num2str(relres)]);
disp(['Τελική Α-νόρμα του Σφάλματος (Τροποποιημένο PCG): ', num2str(errvec(end))]);
