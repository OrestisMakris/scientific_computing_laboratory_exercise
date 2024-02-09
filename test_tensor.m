% Initialization
clear;
tol = 1e-8;
nd = 3;
rng(0);
err = zeros(1, nd + 2);
ndim = [2, 3, 4];
Atemp = randi(5, ndim);
Btemp = randi(4, ndim);
X = randi([-1, 1], max(ndim), 1);
A = tensor(Atemp);
B = tensor(Btemp);

% TTV multiplication check
try
    for k = 1:nd
        err(k) = norm(ttv_myid_1084516(Atemp, X(1:ndim(k), 1), k) - double(ttv(A, X(1:ndim(k), 1), k)));
    end
    assert(max(err) < tol, 'TTM multiplication fails');
catch ME1
    disp(ME1.message);
end

% TTT outer multiplication check
try
    err(nd + 1) = norm(tensor(ttt_myid_1084516(Atemp, Btemp )) - ttt(A, B));
    assert(err(nd + 1) < tol, 'TTT outer multiplication fails');
catch ME2
    disp(ME2.message);
end

% TTT inner product check
try
    err(nd + 2) = abs(ttt_myid_1084516(Atemp, Btemp , 'all') - double(ttt(A, B, [1:nd])));
    assert(err(nd + 2) < tol, 'TTT inner product fails');
catch ME3
    disp(ME3.message);
end

% Display error messages if exceptions occurred
if exist('ME1', 'var')
    disp(ME1.message);
end
if exist('ME2', 'var')
    disp(ME2.message);
end
if exist('ME3', 'var')
    disp(ME3.message);
end

