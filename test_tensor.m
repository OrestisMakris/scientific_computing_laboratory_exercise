clear
nd = 3;
rng(0);
err = zeros(1, nd + 2);
ndim = [2, 3, 4];
Atemp = randi(5, ndim);
X = randi([-1, 1], max(ndim), 1);
A = tensor(Atemp);

for k = 1:nd
    result_ttv_myid = ttv_myid(Atemp, X(1:ndim(k)), k);
    result_ttv_toolbox = double(ttv(A, X(1:ndim(k)), k));
    
    % Display the results
    disp(['Result from ttv_myid (dimension ', num2str(k), '):']);
    disp(result_ttv_myid);
    
    disp(['Result from ttv (Tensor Toolbox, dimension ', num2str(k), '):']);
    disp(result_ttv_toolbox);
    
    % Calculate and display the error
    err(k) = norm(result_ttv_myid - result_ttv_toolbox);
    disp(['Error for dimension ', num2str(k), ': ', num2str(err(k))]);
end

disp('Errors:');
disp(err);
