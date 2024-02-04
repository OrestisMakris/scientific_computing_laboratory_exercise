function Y = ttv_myid(X, V, N)
    % Check if the specified dimension N is valid
    if N < 1 || N > ndims(X)
        error('Invalid dimension specified.');
    end
    
    % Get the size of X
    sz = size(X);
    
    % Check if the size of V is compatible with the specified dimension
    if numel(V) ~= sz(N)
        error('Size mismatch between dimension N and the size of V.');
    end
    
    % Permute the dimensions of X to bring the specified dimension to the front
    perm_order = [N, 1:N-1, N+1:ndims(X)];
    X_permuted = permute(X, perm_order);
    
    % Reshape X to a 2D matrix for multiplication
    X_reshaped = reshape(X_permuted, sz(N), []);
    
    % Perform the multiplication (making sure the dimensions are correct)
    result_2d = X_reshaped.' * V;  % Transpose to match dimensions
    
    % Reshape the result back to a 2D matrix without singleton dimensions
    Y = reshape(result_2d.', size(X, setdiff(1:ndims(X), N)));
end
