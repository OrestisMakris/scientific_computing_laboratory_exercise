% ***************************************************************************************************************
% DESCLAMER THIS CODE WAS DEVELME WITH THE EXTENSIVE HELP OF DEEPSEEK R1 PROVIDING THE CODE I DEVELOPED FOR THE test_tridReduct_1084516.m
% ****************************************************************************************************************
% Στόχος: Εξέταση της συμπεριφοράς των νορμών των off-diagonals (υποδιαιτήτων a και c) 
% κατά τα βήματα της παράλληλης αναγωγής και προσπάθεια εκτίμησης προσέγγισης 
% χωρίς να εκτελεστούν όλα τα O(log₂(n)) βήματα.
%
% Σε αυτό το πείραμα χρησιμοποιούμε το ίδιο τριδιαγώνιο μητρώο που έχουμε και στην 
% δοκιμή, αλλά διακόπτουμε την αναγωγή νωρίτερα (partial reduction). 

%% Παράμετροι πειράματος
n = 1024;                    % μέγεθος του συστήματος (πρέπει να είναι δύναμη του 2)
maxLevels = log2(n);         % πλήθος πλήρων επιπέδων αναγωγής
partialLevel = maxLevels - 2; % διακοπή 2 επίπεδα πριν το τέλος

% Κατασκευή του τριδιαγώνιου μητρώου σε μορφή LAPACK:
A_band = zeros(3, n);
A_band(2,:) = 4;             % κύρια διαγώνιος
A_band(1,2:end) = -1;         % υπερδιαγώνιος
A_band(3,1:end-1) = -1;       % υποδιαγώνιος

% Δημιουργία πλήρους μητρώου για έλεγχο:
A_full = toeplitz([4, -1, zeros(1, n-2)]);

% Κατασκευή δεξιού μέλους:
B = A_full * ones(n, 1);

%% Αποθήκευση καταστάσεων για κάθε επίπεδο
levels = maxLevels;
A_cell = cell(levels+1,1);  % αποθήκευση a
B_cell = cell(levels+1,1);  % αποθήκευση b
C_cell = cell(levels+1,1);  % αποθήκευση c
D_cell = cell(levels+1,1);  % αποθήκευση d

% Αρχική διάσπαση των διαγωνίων
a = [0; A_band(3,1:end-1).'];      % υποδιαγώνιος (με πρώτο στοιχείο 0)
b = A_band(2,:).';                 % κύρια διαγώνιος
c = [A_band(1,2:end).'; 0];         % υπερδιαγώνιος (με τελευταίο στοιχείο 0)
d = B;                            % δεξί μέλος

A_cell{1} = a;
B_cell{1} = b;
C_cell{1} = c;
D_cell{1} = d;

fprintf('Επίπεδο 0: ||a||_inf = %.2e , ||c||_inf = %.2e\n', norm(a, inf), norm(c, inf));

%% Προώθηση (forward reduction) μέχρι το partial level
for lev = 1:maxLevels
    a_prev = A_cell{lev};
    b_prev = B_cell{lev};
    c_prev = C_cell{lev};
    d_prev = D_cell{lev};
    N = length(b_prev);
    
    % Επιλογή even δεικτών
    even = (2:2:N).';
    newN = length(even);
    
    a_new = zeros(newN,1);
    b_new = zeros(newN,1);
    c_new = zeros(newN,1);
    d_new = zeros(newN, size(d_prev,2));
    
    idx_left = even - 1;
    idx_right = even + 1;
    
    hasLeft = (even > 1);
    hasRight = (even < N);
    
    % Υπολογισμός νέας κύριας διαγωνίου
    b_new(:) = b_prev(even);
    if any(hasLeft)
        b_new(hasLeft) = b_new(hasLeft) - a_prev(even(hasLeft)).*( c_prev(idx_left(hasLeft)) ./ b_prev(idx_left(hasLeft)) );
    end
    if any(hasRight)
        b_new(hasRight) = b_new(hasRight) - c_prev(even(hasRight)).*( a_prev(idx_right(hasRight)) ./ b_prev(idx_right(hasRight)) );
    end
    
    % Υπολογισμός νέων off-diagonals
    a_new(hasLeft) = - a_prev(even(hasLeft)).*( a_prev(idx_left(hasLeft)) ./ b_prev(idx_left(hasLeft)) );
    c_new(hasRight) = - c_prev(even(hasRight)).*( c_prev(idx_right(hasRight)) ./ b_prev(idx_right(hasRight)) );
    
    % Δεξί μέλος
    d_new(:,:) = d_prev(even,:);
    if any(hasLeft)
        d_new(hasLeft,:) = d_new(hasLeft,:) - a_prev(even(hasLeft)).*( d_prev(idx_left(hasLeft),:) ./ b_prev(idx_left(hasLeft)) );
    end
    if any(hasRight)
        d_new(hasRight,:) = d_new(hasRight,:) - c_prev(even(hasRight)).*( d_prev(idx_right(hasRight),:) ./ b_prev(idx_right(hasRight)) );
    end
    
    % Αποθήκευση του επιπέδου
    A_cell{lev+1} = a_new;
    B_cell{lev+1} = b_new;
    C_cell{lev+1} = c_new;
    D_cell{lev+1} = d_new;
    
    % Εκτύπωση νορών για το επίπεδο
    fprintf('Επίπεδο %d: ||a||_inf = %.2e , ||c||_inf = %.2e\n', lev, norm(a_new, inf), norm(c_new, inf));
    
    % Διακοπή μετά το partialLevel
    if lev == partialLevel
        fprintf('\nΔιακοπή μετά από επίπεδο %d για partial λύση.\n\n', lev);
        break;
    end
end

%% Back–substitution με το υπάρχον (partial) επίπεδο
currentLevel = lev;  % το τελευταίο επίπεδο που επεξεργαστήκαμε
x_level = D_cell{currentLevel+1} ./ B_cell{currentLevel+1};

for L = currentLevel:-1:1
    a_prev = A_cell{L};
    b_prev = B_cell{L};
    c_prev = C_cell{L};
    d_prev = D_cell{L};
    N = length(b_prev);
    
    x_curr = zeros(N, size(d_prev,2));
    even = (2:2:N).';
    x_curr(even,:) = x_level;
    
    odd = (1:2:N).';
    idx_left = odd - 1;
    idx_right = odd + 1;
    
    x_temp = d_prev(odd,:);
    maskL = odd > 1;
    if any(maskL)
        x_temp(maskL,:) = x_temp(maskL,:) - a_prev(odd(maskL)).*x_level((odd(maskL)-1)/2, :);
    end
    maskR = odd < N;
    if any(maskR)
        x_temp(maskR,:) = x_temp(maskR,:) - c_prev(odd(maskR)).*x_level((odd(maskR)+1)/2, :);
    end
    x_curr(odd,:) = x_temp ./ b_prev(odd);
    x_level = x_curr;
end
x_approx = x_level;

%% Λύση με πλήρη αναγωγή (backslash)
x_full = A_full \ B;

%% Σύγκριση αποτελεσμάτων
diff_partial = norm(x_approx - x_full, inf);
fprintf('Μέγιστη διαφορά μεταξύ της partial λύσης (με %d επίπεδα) και της πλήρους λύσης: %.2e\n', partialLevel, diff_partial);
fprintf('Αν η διαφορά είναι μεγάλη (και όχι το πολύ τάξης O(κ(A)*eps)) τότε ο κώδικάς σας δεν λειτουργεί σωστά.\n');