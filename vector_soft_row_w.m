function Y = vector_soft_row_w(X,tau)

NU = sqrt(sum(X.^2,2));
A = max(0, NU-tau ./ (NU + eps));
Y = repmat((A./(A+tau ./ (NU + eps))),1,size(X,2)).* X;
