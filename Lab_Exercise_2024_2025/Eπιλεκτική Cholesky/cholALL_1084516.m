function [ld, ls] = cholALL_1084516(u,d)
% cholALL_ID: Υπολογίζει τα στοιχεία λ(j,j) και λ(j,j-1) του
%             Cholesky παράγοντα L του A = diag(d) + u*u'
%             για όλα τα j = 1:n.
%
% Έξοδοι:
%   ld : Διανύσμα μεγέθους n με τα στοιχεία της διαγωνίου (λ(j,j))
%   ls : Διανύσμα μεγέθους n με τα στοιχεία της υποδιαγωνίου (λ(j,j-1)),
%        όπου ls(1) = 0.

n = length(u);
ld = zeros(n,1);
ls = zeros(n,1);
L = zeros(n,n);  % Θα αποθηκεύουμε τα απαιτούμενα στοιχεία του L

for i = 1:n
    % Υπολογισμός των μη διαγωνίων στοιχείων της i-οστής γραμμής:
    for k = 1:i-1
        temp = 0;
        for s = 1:k-1
            temp = temp + L(i,s)*L(k,s);
        end
        L(i,k) = ( u(i)*u(k) - temp ) / L(k,k);
    end
    % Υπολογισμός του διαγωνίου στοιχείου:
    temp = 0;
    for s = 1:i-1
        temp = temp + L(i,s)^2;
    end
    L(i,i) = sqrt( d(i) + u(i)^2 - temp );
    ld(i) = L(i,i);
    if i == 1
        ls(i) = 0;
    else
        ls(i) = L(i,i-1);
    end
end
end
