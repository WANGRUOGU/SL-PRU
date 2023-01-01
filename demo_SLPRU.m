%% Loading E.coli reference images and an unknown mixture image
clear variables
load('EcoliImages.mat')
%% Endmembers extraction
E_Ecoli = PNMF(ref_Ecoli);
%% Tuning parameters
lam1 = [0, 10.^(-2:1)];
lam2 = [0, 10.^(-2:1)];
sz_sam = 1e2;
[l1_opt, l2_opt] = Tuning(ref_Ecoli, E_Ecoli, Fluorophores, lam1, lam2, sz_sam);
%% Unmixing a real mixed biological image
A = SLPRU_Full(mixture, E_Ecoli, l1_opt, l2_opt);
%% Estimated abundance matrices
sz_ref = size(ref_Ecoli, 2);
sz_row = floor(sqrt(sz_ref));
sz_col = ceil(sz_ref / sz_row);
figure
for r = 1:sz_ref
    subplot(sz_row, sz_col, r), imshow(A(:, :, r))
    title(Fluorophores{r})
end