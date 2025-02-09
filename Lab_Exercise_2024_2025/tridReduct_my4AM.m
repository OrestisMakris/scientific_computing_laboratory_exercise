function X = tridReduct_my4AM(A, B)
    % tridReduct_my4AM: Λύση τριδιαγώνιων συστημάτων με Κυκλική Αναγωγή (Cyclic Reduction)
    % A: τριδιαγώνιο μητρώο στην ειδική του μορφή [L, D, U]
    % B: πίνακας των δεξιών μελών (κάθε στήλη ένα διαφορετικό σύστημα)
    % X: λύση του συστήματος

    [n, m] = size(B);
    L = A(:, 1); % Υποδιαγώνιος
    D = A(:, 2); % Κύρια διαγώνιος
    U = A(:, 3); % Υπερδιαγώνιος
    
    log_n = log2(n);
    if mod(n, 2^floor(log_n)) ~= 0
        error('Το n πρέπει να είναι δύναμη του 2 για Κυκλική Αναγωγή.');
    end
    
    % Κυκλική Αναγωγή (CR)
    for step = 1:log_n
        stride = 2^(step-1);
        for i = stride+1:2*stride:n
            factor1 = L(i) / D(i-stride);
            factor2 = U(i-stride) / D(i-stride);
            D(i) = D(i) - factor1 * U(i-stride);
            if i+stride <= n
                D(i) = D(i) - factor2 * L(i+stride);
            end
            B(i, :) = B(i, :) - factor1 * B(i-stride, :);
            if i+stride <= n
                B(i, :) = B(i, :) - factor2 * B(i+stride, :);
            end
        end
    end
    
    % Οπίσθια αντικατάσταση
    X = zeros(n, m);
    X(1:2:end, :) = B(1:2:end, :) ./ D(1:2:end);
    for step = log_n:-1:1
        stride = 2^(step-1);
        for i = stride+1:2*stride:n
            X(i, :) = (B(i, :) - L(i) .* X(i-stride, :));
            if i+stride <= n
                X(i, :) = X(i, :) - U(i) .* X(i+stride, :);
            end
            X(i, :) = X(i, :) ./ D(i);
        end
    end
end