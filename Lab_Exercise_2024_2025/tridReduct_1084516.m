function x = tridReduct_1084516(A, B)
% tridReduct_1084516: Λύση τριδιαγώνιου συστήματος με παράλληλη αναγωγή
%
% Είσοδοι:
%   A : Τριδιαγώνιο μητρώο σε μορφή LAPACK (3×n πίνακας)
%   B : Δεξιά μέλη, πίνακας n×m (m≥1)
%
% Έξοδος:
%   x : Λύση του συστήματος, πίνακας n×m
%
% Σημείωση: Θεωρούμε n = 2^k ώστε να εκτελεστεί log₂(n) αναγωγικών βημάτων.
%
% Πολυπλοκότητα: Χρησιμοποιώντας διαδικασίες Hadamard (element–wise πολλαπλασιαστές)
% η διαδικασία απαιτεί O(n·log(n)) πράξεις σε σειριακή εκτέλεση (με παράλληλη υλοποίηση O(log(n)) βήματα).

[nRows, n] = size(A);

% Εξαγωγή διαγωνίων (στήλες ως: [υπερδιαγώνιος; κύρια διαγώνιος; υποδιαγώνιος])
a = [0; A(3,1:end-1).'];    % υποδιαγώνιος (στήλη) με πρώτο στοιχείο 0
b = A(2,:).';               % κύρια διαγώνιος (στήλη)
c = [A(1,2:end).'; 0];       % υπερδιαγώνιος (στήλη) με τελευταίο στοιχείο 0

% Δεξί μέλος
d = B;  % n×m

% Αποθήκευση συντελεστών για back–substitution
levels = log2(n);
A_cell = cell(levels+1,1);  % αποθήκευση υποδιαγωνίων
B_cell = cell(levels+1,1);  % αποθήκευση κύριων διαγωνίων
C_cell = cell(levels+1,1);  % αποθήκευση υπερδιαγωνίων
D_cell = cell(levels+1,1);  % αποθήκευση δεξιών μελών

A_cell{1} = a;
B_cell{1} = b;
C_cell{1} = c;
D_cell{1} = d;

% ****************************
% ΠΡΟΣΩΡΙΝΗ ΟΜΟΓΕΝΟΠΟΙΗΣΗ (FORWARD REDUCTION)
% ****************************
% Σε κάθε βήμα μειώνεται το μέγεθος του συστήματος κατά 2,
% χρησιμοποιώντας vectorized ενημερώσεις τύπων:
%
%   b_even_new = b_even - a_even .* (c_odd ./ b_odd) - c_even .* (a_odd2 ./ b_odd2)
%   a_even_new = - a_even .* (a_odd ./ b_odd)
%   c_even_new = - c_even .* (c_odd2 ./ b_odd2)
%   d_even_new = d_even - a_even .* (d_odd ./ b_odd) - c_even .* (d_odd2 ./ b_odd2)
%
% όπου για κάθε even-index i στο τρέχον επίπεδο:
%   - o₁ = i-1 αν υπάρχει (για ορισμό: a(i) συζεύγεται με τον αριστερό γείτονα)
%   - o₂ = i+1 αν υπάρχει (για τον δεξιό γείτονα)
%
% Οι ενημερώσεις γίνονται με element–wise (Hadamard) πράξεις.
for lev = 1:levels
    a_prev = A_cell{lev};
    b_prev = B_cell{lev};
    c_prev = C_cell{lev};
    d_prev = D_cell{lev};
    N = length(b_prev);
    
    % Επιλέγουμε τους "even" δείκτες στο τρέχον σύστημα
    even = (2:2:N).';
    newN = length(even);
    
    % Προετοιμασία νέων διανυσμάτων (για τους even δείκτες)
    a_new = zeros(newN,1);
    b_new = zeros(newN,1);
    c_new = zeros(newN,1);
    d_new = zeros(newN, size(d_prev,2));
    
    % Για κάθε even-index i στο επίπεδο, οι γείτονες:
    % ο αριστερός δείκτης: iL = i-1  (αν υπάρχει)
    % ο δεξιός δείκτης:   iR = i+1  (αν υπάρχει)
    idx_left = even - 1;  
    idx_right = even + 1;
    
    % Εφαρμόζουμε element–wise τους τύπους όπου υπάρχουν και οι δύο γείτονες:
    % Δημιουργούμε μάσκες για τις περιπτώσεις οριακών τιμών.
    hasLeft = (even > 1);
    hasRight = (even < N);
    
    % b_new: αρχικά ίσο με b(even)
    b_new(:) = b_prev(even);
    % Αφαιρούμε τους όρους μόνο όπου υπάρχουν οι γείτονες
    b_new(hasLeft) = b_new(hasLeft) - a_prev(even(hasLeft)) .* ( c_prev(idx_left(hasLeft)) ./ b_prev(idx_left(hasLeft)) );
    b_new(hasRight) = b_new(hasRight) - c_prev(even(hasRight)) .* ( a_prev(idx_right(hasRight)) ./ b_prev(idx_right(hasRight)) );
    
    % a_new: μόνο αν υπάρχει αριστερός γείτονας
    a_new(hasLeft) = - a_prev(even(hasLeft)) .* ( a_prev(idx_left(hasLeft)) ./ b_prev(idx_left(hasLeft)) );
    % c_new: μόνο αν υπάρχει δεξιός γείτονας
    c_new(hasRight) = - c_prev(even(hasRight)) .* ( c_prev(idx_right(hasRight)) ./ b_prev(idx_right(hasRight)) );
    
    % d_new: ενημέρωση δεξιού μέλους (για όλες τις στήλες)
    d_new(:,:) = d_prev(even,:);
    if any(hasLeft)
        d_new(hasLeft,:) = d_new(hasLeft,:) - a_prev(even(hasLeft)) .* ( d_prev(idx_left(hasLeft),:) ./ b_prev(idx_left(hasLeft)) );
    end
    if any(hasRight)
        d_new(hasRight,:) = d_new(hasRight,:) - c_prev(even(hasRight)) .* ( d_prev(idx_right(hasRight),:) ./ b_prev(idx_right(hasRight)) );
    end
    
    % Αποθήκευση του νέου (μειωμένου) συστήματος στο επόμενο επίπεδο
    A_cell{lev+1} = a_new;
    B_cell{lev+1} = b_new;
    C_cell{lev+1} = c_new;
    D_cell{lev+1} = d_new;
end

% Στο τελευταίο επίπεδο, το σύστημα είναι διαγώνιο:
x_reduced = D_cell{end} ./ B_cell{end};

% ****************************
% BACK–SUBSTITUTION
% ****************************
% Επαναφορά των λύσεων στα πλήρη συστήματα του κάθε επιπέδου.
% Σε κάθε βήμα, για τους "odd" δείκτες, χρησιμοποιούμε τη σχέση:
%
%   x(i) = ( d(i) - a(i)*x(i-1) - c(i)*x(i+1) ) ./ b(i)
%
% οι ενημερώσεις γίνονται και εδώ με element–wise πράξεις.
% BACK–SUBSTITUTION
x_level = x_reduced;
for lev = levels:-1:1
    a_prev = A_cell{lev};
    b_prev = B_cell{lev};
    c_prev = C_cell{lev};
    d_prev = D_cell{lev};
    N = length(b_prev);
    
    % Δημιουργούμε το διάνυσμα λύσεων για το τρέχον επίπεδο
    x_curr = zeros(N, size(d_prev,2));
    
    % Οι λύσεις για τους even δείκτες έχουν ήδη υπολογιστεί στο επίπεδο
    even = (2:2:N).';
    x_curr(even,:) = x_level;
    
    % Υπολογισμός για τους odd δείκτες
    odd = (1:2:N).';
    idx_left = odd - 1;  
    idx_right = odd + 1;
    
    % Αρχικά, δίνει την τιμή d_prev(odd)
    x_temp = d_prev(odd,:);
    % Απουσία αριστερού γείτονα για το πρώτο στοιχείο
    maskL = odd > 1;
    if any(maskL)
        x_temp(maskL,:) = x_temp(maskL,:) - a_prev(odd(maskL)) .* x_level((odd(maskL)-1)/2, :);
    end
    % Για δεξί γείτονα, χρησιμοποιούμε τη διωνυμική αντιστοίχιση (i+1)/2
    maskR = odd < N;
    if any(maskR)
        x_temp(maskR,:) = x_temp(maskR,:) - c_prev(odd(maskR)) .* x_level((odd(maskR)+1)/2, :);
    end
    x_curr(odd,:) = x_temp ./ b_prev(odd);
    
    x_level = x_curr;
end
x = x_level;
end