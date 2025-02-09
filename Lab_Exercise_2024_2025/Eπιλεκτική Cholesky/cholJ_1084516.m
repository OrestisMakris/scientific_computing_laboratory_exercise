function [ljj, ljj_1] = cholJ_ID(u,d,j)
% cholJ_ID: Υπολογίζει τα στοιχεία λ(j,j) και λ(j,j-1) του Cholesky παράγοντα
%           του A = diag(d) + u*u' για τον υποπίνακα 1:j,1:j.
%
% Είσοδοι:
%   u : Διάνυσμα στο ℝ^n
%   d : Διάνυσμα στο ℝ^n με d(i) > 0 για κάθε i
%   j : Δείκτης (1 ≤ j ≤ n)
%
% Έξοδοι:
%   ljj   : Το στοιχείο L(j,j)
%   ljj_1 : Το στοιχείο L(j,j-1) (για j = 1, επιστρέφεται 0)

n = length(u);
if j < 1 || j > n
    error('Ο δείκτης j πρέπει να 1 ≤ j ≤ n');
end

% Αρχικοποίηση του πίνακα L για τους πρώτους j γραμμές και στήλες.
L = zeros(j,j);

for i = 1:j
    % Υπολογισμός των μη διαγωνίων στοιχείων της i-οστής γραμμής:
    for k = 1:i-1
        temp = 0;
        for s = 1:k-1
            temp = temp + L(i,s)*L(k,s);
        end
        % Για i ≠ k έχουμε A(i,k) = u(i)*u(k)
        L(i,k) = ( u(i)*u(k) - temp ) / L(k,k);
    end
    % Υπολογισμός του διαγωνίου στοιχείου:
    temp = 0;
    for s = 1:i-1
        temp = temp + L(i,s)^2;
    end
    L(i,i) = sqrt( d(i) + u(i)^2 - temp );
end

ljj = L(j,j);
if j == 1
    ljj_1 = 0;
else
    ljj_1 = L(j,j-1);
end
end
