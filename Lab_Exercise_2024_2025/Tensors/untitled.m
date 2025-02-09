% Σε αυτήν την έκδοση, χρησιμοποιούμε εικόνες που είναι ήδη αποθηκευμένες τοπικά.
% Εδώ θα χρησιμοποιήσουμε τις εικόνες "egkefalos.jpg" και "p5boat.jpg".

local_filenames = {'egkefalos.jpg', 'p5boat.jpg'};

% Αρχικοποίηση τυχαίας γεννήτριας αριθμών
rng(42);

% Επεξεργασία κάθε τοπικά αποθηκευμένης εικόνας
for i = 1:length(local_filenames)
    % Διαβάστε την εικόνα από το αρχείο.
    % Δεν μετατρέπουμε την εικόνα σε grayscale ώστε να διατηρηθούν τα χρώματα.
    try
        img = imread(local_filenames{i});
    catch ME
        fprintf('Error reading image %s: %s\n', local_filenames{i}, ME.message);
        continue;
    end
    
    % Δημιουργία τανυστή από την εικόνα.
    % Για έγχρωμες jpg, το img είναι ένας τρισδιάστατος πίνακας (rows x cols x 3).
    % Μετατρέπουμε τον πίνακα σε αντικείμενο tensor (Tensor Toolbox) για να αποφευχθεί το σφάλμα στον υπολογισμό του norm μέσα στην cp_als.
    T = tensor(double(img));
    original_image = double(img);
    
    ranks = [];
    f_deltas = [];
    approximations = {};
    
    % Υπολογίστε τις προσεγγίσεις για τα διάφορα f-delta thresholds
    for f_delta_threshold = [0.5, 0.7, 0.9]
        k = 1;
        while true
            % Εκτέλεση της διάσπασης CP χρησιμοποιώντας ALS για τον τρέχοντα βαθμό (k)
            % (cp_als αναμένει ένα tensor αντικείμενο, οπότε χρησιμοποιούμε το T.)
            P = cp_als(T, k, 'printitn', 0);
            % Ανακατασκευή του τανυστή από την διάσπαση
            reconstructed_tensor = double(full(P));
            % Υπολογισμός του f-delta ως λόγος των νόρμων
            f_delta = norm(original_image(:) - reconstructed_tensor(:)) / norm(original_image(:));
            
            if f_delta > f_delta_threshold
                ranks(end + 1) = k;
                f_deltas(end + 1) = f_delta;
                approximations{end + 1} = reconstructed_tensor;
                break;
            end
            k = k + 1;
        end
    end
    
    % Οπτικοποίηση της αρχικής εικόνας και των προσεγγίσεων
    figure;
    subplot(2, 2, 1);
    imshow(uint8(original_image));
    title('Original Image');
    
    for j = 1:3
        subplot(2, 2, j + 1);
        imshow(uint8(approximations{j}));
        title(sprintf('k = %d, f-delta = %.2f', ranks(j), f_deltas(j)));
    end
    
    % Υπολογισμός και εμφάνιση του κόστους αποθήκευσης για κάθε αναπαράσταση
    for j = 1:length(ranks)
        k = ranks(j);
        % Επανεκτέλεση της διάσπασης CP για τον συγκεκριμένο βαθμό (k)
        P = cp_als(T, k, 'printitn', 0);
        % Ανάκτηση της αναπαράστασης Kruskal (cell array με τους παράγοντες)
        kruskal_representation = P.U;
        original_storage_cost = numel(original_image);
        kruskal_storage_cost = sum(cellfun(@numel, kruskal_representation));
        full_array_storage_cost = numel(original_image);
        
        fprintf('Για k = %d και f-delta = %.2f:\n', k, f_deltas(j));
        fprintf('Κόστος αποθήκευσης (πλήρης τανυστής): %d\n', original_storage_cost);
        fprintf('Κόστος αποθήκευσης (αναπαράσταση Kruskal): %d\n', kruskal_storage_cost);
        fprintf('Κόστος αποθήκευσης (πολυδιάστατο array): %d\n\n', full_array_storage_cost);
    end
end