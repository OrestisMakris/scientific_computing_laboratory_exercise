% Τανυστική Ανάλυση Εικόνων με CP-ALS χρησιμοποιώντας το Tensor Toolbox

% Αρχικοποίηση τυχαίων αριθμών

rng(1234);

% Φόρτωση εικόνων και μείωση μεγέθους για ταχύτερη επεξεργασία xrish
% toolbox gia image procesing 

img1 = imresize(imread('egkefalos.jpg'), 0.2);

img2 = imresize(imread('p5boat.jpg'), 0.2);



% Μετατροπή εικόνων σε τανυστές

tensor1 = tensor(double(img1) / 255);
 

tensor2 = tensor(double(img2) / 255);

% Ορισμός κατωφλίων για f-delta

f_thresholds = [0.5,0.7,0.9];

k_values = zeros(length(f_thresholds), 2);

approx_tensors = cell(length(f_thresholds), 2); % To store the approximated tensors

% Συνάρτηση για εύρεση του k για κάθε f-delta
for img_idx = 1:2

    tensor_var = eval(['tensor',  num2str(img_idx ) ]);  % Dynamicaly select tensor  1 or tensor 2

    dims = size(tensor_var);

    for t = 1:length(f_thresholds)

        k = 1;

        f_delta = 0;


        while f_delta < f_thresholds(t)
            
            % Διάσπαση CP-ALS toolbox edo 

            M = cp_als(tensor_var, k, 'printitn', 0, 'tol', 1e-5, 'maxiters', 80);  %parametru gia tolmaxiters
            approx_tensorr = full(M);
            
            % Υπολογισμός f-delta χρησιμοποιώντας την Frobenius norm
            f_delta = 1 - norm(double(tensor_var) - double(approx_tensorr), 'fro') / norm(double(tensor_var), 'fro');
            
            if f_delta >= f_thresholds(t)
                k_values(t, img_idx) = k;
                approx_tensors{t, img_idx} =  approx_tensor;                 
                break;
            end
            k = k + 1;
            
            % Προσθήκη ελέγχου για αποφυγή απεριόριστης εκτέλεσης

            if k > 400  % Adjusted maximum k value

                warning('Error reached maximum k value without meeting threshold for image %d at threshold  %.1f', img_idx, f_thresholds(t));

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

    fprintf('Αποθηκευση πλιρους τανυστη σε: %d στοιχεία\n', original_storage);
    
    for t = 1:length(f_thresholds)

        k = k_values(t, img_idx);

        storage_kruskal = k * sum(dims);

        fprintf('Αποθηκευση Kruscal (k=%d): %d στοιχεiα\n', k, storage_kruskal);
    end
    fprintf('\n');
end


figure;

for img_idx = 1:2

    tensor_var = eval(['tensor', num2str(img_idx)]);  % Dynamically select tensor1 or tensor2
    subplot(4, 2, (img_idx - 1) * 4 + 1);
    imshow(double(tensor_var));
    title('Αρχική Εικόνα');
    
    for t = 1:length(f_thresholds)

        k = k_values(t, img_idx);
        approx_tensorr = approx_tensors{t, img_idx}; 
        subplot(4, 2, (img_idx - 1) * 4 + 1 + t);
        imshow(double(approx_tensorr));
        title(sprintf('k=%d, f-delta=%.1f', k, f_thresholds(t)));
    end
end
