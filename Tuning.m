function [l1_PR, l2_PR, l1_NLS, l2_NLS, l_SNLS, prop]...
    = Tuning(ref, E, E_name, lam1, lam2, sz_sam)

sz_l1 = length(lam1);
sz_l2 = length(lam2);
[C, sz_ref] = size(E);
prop_NLS_temp = zeros(sz_l1, sz_l2, sz_ref, sz_sam);
prop_PR_temp = zeros(sz_l1, sz_l2, sz_ref, sz_sam);
for r = 1:sz_ref
    M = ref{r};
    ind = find(any(M, 2));
    sz = sqrt(size(M, 1));
    T = reshape(M, sz, sz, C);
    [row, col] = ind2sub([sz sz], ind);
    ind = find(row > 1 & row < sz & col > 1 & col < sz);
    row = row(ind);
    col = col(ind);
    ind_sam = randsample(length(ind), sz_sam);
    for i = 1:sz_sam
        ind_x = (row(ind_sam(i))-1):(row(ind_sam(i))+1);
        ind_y = (col(ind_sam(i))-1):(col(ind_sam(i))+1);
        M = reshape(T(ind_x, ind_y, :), length(ind_x)*length(ind_y), C);
        for l1 = 1:sz_l1
            for l2 = 1:sz_l2
                [~, A] = SLPRU(M.', E, lam1(l1), lam2(l2));
                prop_PR_temp(l1, l2, r, i) = A(r) / sum(A);
                [~, A] = SLNLS(M.', E, lam1(l1), lam2(l2));
                prop_NLS_temp(l1, l2, r, i) = A(r) / sum(A);
            end
        end
        fprintf('Pixels of Ref %s unmixed: %2d/%2d \n', E_name{r}, i, sz_sam);
    end
end
% optimal parameters for SL-PRU 
prop = mean(prop_PR_temp, 4, 'omitnan');
prop_min = min(prop, [], 3);
[~, ind] = max(prop_min, [], 'all');
[ind1, ind2] = ind2sub(size(prop_min), ind);
l1_PR = lam1(ind1);
l2_PR = lam1(ind2);
% optimal parameters for S-NLS
prop = mean(prop_NLS_temp, 4, 'omitnan');
prop_SNLS = reshape(min(prop(1, :, :), [], 3), sz_l2, 1);
[~, ind_SNLS] = max(prop_SNLS);
l_SNLS = lam2(ind_SNLS);
% optimal parameters for SL-NLS
prop_min = min(prop, [], 3);
[~, ind] = max(prop_min, [], 'all');
[ind1, ind2] = ind2sub(size(prop_min), ind);
l1_NLS = lam1(ind1);
l2_NLS = lam2(ind2);