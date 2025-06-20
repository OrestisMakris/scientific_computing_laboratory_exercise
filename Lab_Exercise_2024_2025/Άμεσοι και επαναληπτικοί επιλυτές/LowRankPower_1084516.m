function [U,V] = LowRankPower_id(A, u, v, k)
% LowRankPower_id
% Υπολογίζει τους παράγοντες U, V (με διαστάσεις n×k) για τους οποίους
% ισχύει: (A+u*v')^k = A^k + U*V'.
%
% Εισόδους:
%   A  : τριδιαγώνιος, αραιός και ΣΘΟ πίνακας (n×n)
%   u,v: στήλες (n×1)
%   k  : θετικός ακέραιος, με 1 < k < n
%
% Εξόδους
% :
%   U,V: πίνακες διαστάσεων n×k, τέτοιοι ώστε:
%         (A+u*v')k = A^ k + U*V'
%
% Η ιδέα είναι να χρησιμοποιήσουμε την ταυτότητα:
%   (A+u*v')^k - A^k = sum_{j= 0}^{k-1} (A+u*v")^j * u * (A^ (k-1-j]*v)^ T.
%
% Υπολογισμοί γίνονται με επαναληπτικό πολλαπλασιασμό  (matrix–free)
% ώστε να εκμεταλλευτούμε την αραιότητα του A και τ η δομή της rank–1 ενημέρωσης.

n = length(u);
U = zeros(n, k);
V = zeros(n, k);


U(:,1) = u;  % j = 0

for j = 1:k-1

    % Ο πολλαπλασιασμός με (A+u*v') γίνεται ως:
    %   y = A*x + u*(v'*x)

    U(:,j+1) = A * U(:,j) + u * (v' * U(:,j));
end


% Υπολογισμός V: V(:,j+1) = A^(k-1-j)*v, για j = 0,...,k-1.
% Για j = k-1 έχουμε: A^0*v = v.

V(:,k) = v;
for j = k-1:-1:1
    V(:,j) = A * V(:,j+1);
end


end
