% Τανυστική Ανάλυση Εικόνων με CP-ALS χρησιμοποιώντας το Tensor Toolbox

% Αρχικοποίηση τυχαίων αριθμών
rng(1234);

% Φόρτωση εικόνων
img1 = imread('egkefalos.jpg');
img2 = imread('p5boat.jpg');

% Μετατροπή εικόνων σε τανυστές
tensor1 = tensor(double(img1) / 255);
tensor2 = tensor(double(img2) / 255);

% Ορισμός κατωφλίων για f-delta
thresholds = [0.9];
k_values = zeros(length(thresholds), 2);

% Συνάρτηση για εύρεση του k για κάθε f-delta
for img_idx = 1:2
    tensor_var = eval(['tensor', num2str(img_idx)]);  % Dynamically select tensor1 or tensor2
    dims = size(tensor_var);
    for t = 1:length(thresholds)
        k = 1;
        f_delta = 0;
        while f_delta < thresholds(t)
            % Διάσπαση CP-ALS
            M = cp_als(tensor_var, k, 'printitn', 0, 'tol', 1e-4, 'maxiters', 50);
            approx_tensor = full(M);
            
            % Υπολογισμός f-delta χρησιμοποιώντας την Frobenius norm
            f_delta = 1 - norm(double(tensor_var) - double(approx_tensor), 'fro') / norm(double(tensor_var), 'fro');
            
            if f_delta >= thresholds(t)
                k_values(t, img_idx) = k;
                break;
            end
            k = k + 1;
            
            % Προσθήκη ελέγχου για αποφυγή απεριόριστης εκτέλεσης
            if k > 100
                warning('Reached maximum k value without meeting threshold for image %d at threshold %.1f', img_idx, thresholds(t));
                break;
            end
        end
    end
end

% Σύγκριση κόστους αποθήκευσης
for img_idx = 1:2
    tensor_var = eval(['tensor', num2str(img_idx)]);  % Dynamically select tensor1 or tensor2
    original_storage = numel(tensor_var);
    dims = size(tensor_var);
    
    fprintf('Εικόνα %d:\n', img_idx);
    fprintf('Αποθήκευση πλήρους τανυστή: %d στοιχεία\n', original_storage);
    
    for t = 1:length(thresholds)
        k = k_values(t, img_idx);
        storage_kruskal = k * sum(dims);
        fprintf('Αποθήκευση Kruskal (k=%d): %d στοιχεία\n', k, storage_kruskal);
    end
    fprintf('\n');
end

% Οπτικοποίηση αποτελεσμάτων
figure;
for img_idx = 1:2
    tensor_var = eval(['tensor', num2str(img_idx)]);  % Dynamically select tensor1 or tensor2
    subplot(4, 2, (img_idx - 1) * 4 + 1);
    imshow(double(tensor_var));
    title('Αρχική Εικόνα');
    
    for t = 1:length(thresholds)
        k = k_values(t, img_idx);
        M = cp_als(tensor_var, k, 'printitn', 0, 'tol', 1e-4, 'maxiters', 50);  % Suppress iteration output
        approx_tensor = full(M);
        subplot(4, 2, (img_idx - 1) * 4 + 1 + t);
        imshow(double(approx_tensor));
        title(sprintf('k=%d, f-delta=%.1f', k, thresholds(t)));
    end
end

% Αποθήκευση ως MATLAB Live Script
matlab.internal.liveeditor.openAndConvert('cp_als_analysis.mlx');