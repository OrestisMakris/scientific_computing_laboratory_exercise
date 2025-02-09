function x = LowRankPowerSLVpcg_1084516(A, u, v, k, b)
% LowRankPowerSLVpcg_id
% Όπως η LowRankPowerSLV1_id, αλλά οι λύσεις με το A δίνονται με την pcg.
tol = 1e-9;
maxit = min(size(A,1), 100);
% Προϋπολογισμός: λύση A*x = u μέσω pcg
[x_pc, flag] = pcg(A, u, tol, maxit);
if flag ~= 0
    warning('pcg δεν συγκλίνει για την επίλυση A\\u.');
end
A_inv_u = x_pc;
denom = 1 + v' * A_inv_u;
y = b;
for i = 1:k
    [z, flag] = pcg(A, y, tol, maxit);
    if flag ~= 0
        warning('pcg δεν συγκλίνει για A\\y στην επανάληψη %d.', i);
    end
    correction = (v' * z) / denom;
    y = z - A_inv_u * correction;
end
x = y;
end
