function A = UnmixSpIm(T, E, method, l1, l2)

[C, sz_ref] = size(E);
if ismatrix(T)
    sz_im = sqrt(size(T, 1));
    T = reshape(T, [sz_im, sz_im, C]);
else
    T = im2double(T);
    sz_im = size(T, 1);
end

switch method
    case 'NLS'
        f = @(z) SLNLS(z.', E, 0, 0);
    case 'SNLS'
        f = @(z) SLNLS(z.', E, 0, l2);
    case 'SLNLS'
        f = @(z) SLNLS(z.', E, l1, l2);
    case 'SLPRU'
        f = @(z) SLPRU(z.', E, l1, l2);
end

A = zeros(sz_im, sz_im, sz_ref);
for i = 2:(sz_im-1)
    for j = 2:(sz_im-1)
        if any(T(i, j, :))
            ind_x = max(i-1, 1):min(i+1, sz_im);
            ind_y = max(j-1, 1):min(j+1, sz_im);
            M = reshape(T(ind_x, ind_y, :), length(ind_x)*length(ind_y), C);
            [~, A(i, j, :)] = f(M);
            n = (sz_im-2) * (i-2) + j-1; 
            fprintf('Number of pixels unmixed: %d\n', n);
        end
    end
end
A = A ./ max(A(:));






