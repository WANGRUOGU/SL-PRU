function E = PNMF(ref)

iter = 1e3;
sz_ref = length(ref);
fit = 1;
C = size(ref{1}, 2);
E = rand(C, sz_ref);
for r = 1:sz_ref
    M = ref{r};
    M = M(any(M, 2), :);
    A = rand(size(M, 1), 1);
    for i = 1:iter
        fit_old = fit;
        A = A .* (M*E(:, r)) ./ (A*E(:, r).'*E(:, r));
        E(:, r) = E(:, r) .* (M.'*A) ./ (E(:, r)*(A.')*A);
        E(isnan(E)) = 0;
        E(E==Inf) = 0;
        E(:, r) = E(:, r) ./ max(E(:, r));
        A = A .* max(E(:, r));
        fit = norm(M - A * E(:, r).', 'fro');
        change = abs(fit_old - fit)/fit_old;
        fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', i, fit, change);
        if (change < 1e-8) || (isnan(fit))
            break;
        end
    end
end