function A = SLPRU_Full(T, E, l1, l2)

[C, sz_ref] = size(E);
if ismatrix(T)
    sz_im = sqrt(size(T, 1));
    T = reshape(T, [sz_im, sz_im, C]);
else
    T = im2double(T);
    sz_im = size(T, 1);
end
A = zeros(sz_im, sz_im, sz_ref);
for i = 1:sz_im
    for j = 1:sz_im
        if any(T(i, j, :))
            if l1 == 0
                M = reshape(T(i, j, :), 1, C);
            else
                ind_x = max(i-1, 1):min(i+1, sz_im);
                ind_y = max(j-1, 1):min(j+1, sz_im);
                M = reshape(T(ind_x, ind_y, :), length(ind_x)*length(ind_y), C);
            end
            [A(i, j, :), ~] = SLPR(M.', E, l1, l2);
            fprintf(' Number of pixels unmixed: %2d/%2d lam1 = %.3g lam2 = %.3g\n', ...
                (i-1)*sz_im+j, sz_im*sz_im, l1, l2);
        end
    end
end
A = A ./ max(A(:));