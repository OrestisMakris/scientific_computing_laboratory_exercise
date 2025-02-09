function x = LowRankPowerSLV1_1084516(A, u, v, k, b)
% LowRankPowerSLV1_id
% Επιλύει το σύστημα (A+u*v')^k x = b ως εξής:
%   y_1 = b, και για i = 1,...,k λύνουμε (A+u*v') y_{i+1} = y_i.
%
% Χρησιμοποιεί την Sherman–Morrison φόρμουλα για την επίλυση
% του συστήματος με A+u*v'.
%
% Σημείωση: Λύσεις με το A δίνονται μέσω του backslash, εκμεταλλευόμενοι
% τη σπανιότητα του A.
%
% (A+u*v') x = b  =>  x = A^{-1}b - (A^{-1}u*(v'*A^{-1}b))/(1+v'*A^{-1}u)

% Προϋπολογισμός A^{-1}u
A_inv_u = A \ u;
denom = 1 + v' * A_inv_u;
y = b;
for i = 1:k
    z = A \ y;              % Λύση A*z = y
    correction = (v' * z) / denom;
    y = z - A_inv_u * correction;
end
x = y;
end
