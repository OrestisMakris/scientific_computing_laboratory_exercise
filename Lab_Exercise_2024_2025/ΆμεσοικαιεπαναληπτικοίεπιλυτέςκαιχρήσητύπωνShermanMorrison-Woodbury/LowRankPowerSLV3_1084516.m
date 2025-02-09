function x = LowRankPowerSLV3_1084516(A, u, v, k, b)
% LowRankPowerSLV3_id
% Επιλύει το σύστημα (A+u*v')^k x = b χρησιμοποιώντας:
%   1) Τον υπολογισμό των παραγόντων U, V από τη LowRankPower_id.
%   2) Τον τύπο Sherman–Morrison–Woodbury για το A^k + U*V'.
%
% Για την επίλυση με το A^k χρησιμοποιούμε την ιδέα ότι:
%   A^{-k}b = (A^{-1})^k b,
% δηλαδή λύνοντας επανειλημμένα k συστήματα με το A.

[U, V] = LowRankPower_1084516(A, u, v, k);

% Λύνουμε A^k z = b: εφαρμόζουμε k φορές A^{-1}.
z = b;
for j = 1:k
    z = A \ z;
end

% Υπολογισμός W = A^{-k}*U: για κάθε στήλη U(:,j), λύνεται το σύστημα A^k w = U(:,j).
[~, kcols] = size(U);
W = zeros(size(U));
for j = 1:kcols
    w = U(:,j);
    for i = 1:k
        w = A \ w;
    end
    W(:,j) = w;
end

% Μικρό σύστημα: M = I + V'*W
M = eye(k) + V' * W;
y = M \ (V' * z);
x = z - W * y;
end
