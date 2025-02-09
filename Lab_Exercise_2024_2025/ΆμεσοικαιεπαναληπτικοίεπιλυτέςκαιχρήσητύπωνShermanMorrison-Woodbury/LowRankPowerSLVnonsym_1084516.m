function x = LowRankPowerSLVnonsym_1084516(A, u, v, k, b)
% LowRankPowerSLVnonsym_id
% Επιλύει το σύστημα G x = b, όπου G = (A+u*v')^k, χρησιμοποιώντας
% επαναληπτική μέθοδο τύπου Krylov (GMRES) σε matrix–free τρόπο.
tol = 1e-6;
maxit = 100;
% Ορίζουμε τη συνάρτηση για τον πολλαπλασιασμό: 
Ghandle = @(x) power_mult(x, A, u, v, k);
[x, flag, relres] = gmres(Ghandle, b, [], tol, maxit);
if flag ~= 0
    warning('Η GMRES δεν συγκλίνει: flag = %d, relres = %e', flag, relres);
end
end

function y = power_mult(x, A, u, v, k)
% Υπολογίζει y = (A+u*v')^k * x χωρίς να σχηματίζεται ρητά ο G.
y = x;
for i = 1:k
    y = A * y + u * (v' * y);
end
end
