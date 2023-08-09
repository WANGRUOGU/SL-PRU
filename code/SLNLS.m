function [U, A] = SLNLS(Y, M, lam1, lam2)
norm_y = sqrt(mean(Y(:).^2));
Y = Y/norm_y;
[C, N] = size(Y);
R = size(M, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abundance matrix
U = (M'*M)^(-1)*M'*Y;
U(U < 0) = 0;
% Auxiliary matrices
V1 = M*U;
V2 = U;
V3 = U;
V4 = U;
% Scaled Lagrange multipliers
D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
D4 = V4*0;
% Augmented Lagrange regularization parameter
mu = 1e-2;
% error tolerance
epsilon = 1e-4;
tol = sqrt((3*R+C)*N)*epsilon;
% Maximum number of iterations
iter = 1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequences of subproblems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:iter
    if mod(i, 10) == 1
        V1_old = V1;
        V2_old = V2;
        V3_old = V3;
        V4_old = V4;
    end
    U = ((M'*M+3*eye(R))^-1)*(M'*(V1+D1)+(V2+D2)+(V3+D3)+(V4+D4));
    V1 = 1/(1+mu)*(Y + mu*(M*U - D1));
    [u, s, v] = svd(U - D2,'econ');
    ds = diag(s);
    V2 = u*diag(max(ds-(lam1/mu)*(1./(ds+eps)),ds*0))*v';
    V3 = vector_soft_row_w(U-D3, lam2/mu);
    V4 = max(U - D4, 0);
    D1 = D1 - M*U + V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
    D4 = D4 - U + V4;
    if mod(i, 10) == 1
        res_p = norm([M*U; U; U; U] - [V1; V2; V3; V4], 'fro');
        res_d = mu*norm(M.'*(V1-V1_old)+(V2-V2_old)+(V3-V3_old)...
            +(V4-V4_old),'fro');
        if res_p > 10*res_d
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            D4 = D4/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            D4 = D4*2;
        end
    end
    if (res_p < tol) && (res_d < tol)
        break;
    end
end

U = U * norm_y;
U(U < 0) = 0;
A = U(:,round(N/2));
end
