ns = [100, 500, 1000, 2000];
time_cholALL = zeros(size(ns));
time_chol_builtin = zeros(size(ns));

for idx = 1:length(ns)
    n = ns(idx);
    rng('default');
    d = rand(n,1);
    u = randn(n,1);
    A = diag(d) + u*u';
     

      
    % Χρονισμός της cholALL_ID:
    tic;
    [ld, ls] = cholALL_1084516(u,d);
    time_cholALL(idx) = toc;
    
    % Χρονισμός της ενσωματωμένης chol (υπολογισμός του κάτω τριγωνικού L):
    tic;


    L_builtin = chol(A, 'lower');

    time_chol_builtin(idx) = toc;

end

fprintf('n\tcholALL_ID (sec)\tchol (sec)\n');

for idx = 1:length(ns)
    
    fprintf('%d\t%f\t%f\n', ns(idx), time_cholALL(idx), time_chol_builtin(idx));
end
