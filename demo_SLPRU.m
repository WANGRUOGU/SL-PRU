%% Loading E.coli reference images and an unknown mixture image
clear variables
load EcoliEndMat.mat
%% Parameters
lam_seq = @(lam_ub, lam_lb, lam_range) ...
    [0,exp(log(lam_lb):((log(lam_ub)-log(lam_lb))/lam_range):log(lam_ub))];
lam = lam_seq(10, .1, 8);
sz_sam = 1e0;
%% Endmembers extraction
E_Ecoli = PNMF(ref_Ecoli);
%% Tuning parameters
[l1_PR, l2_PR, l1_NLS, l2_NLS, l_SNLS, prop]...
    = Tuning(ref_Ecoli, E_Ecoli, Fluorophores, lam, lam, sz_sam);
%% Unmixing a real mixed biological image
A_NLS = UnmixSpIm(mixture, E_Ecoli, 'NLS', 0, 0);
A_SNLS = UnmixSpIm(mixture, E_Ecoli, 'SNLS', 0, l_SNLS);
A_SLNLS = UnmixSpIm(mixture, E_Ecoli, 'SLNLS', l1_NLS, l2_NLS);
A_SLPRU = UnmixSpIm(mixture, E_Ecoli, 'SLPRU', l1_PR, l2_PR);
%% Estimated abundance matrices
FigSpIm(mixture, A_NLS, A_SNLS, A_SLNLS, A_SLPRU, Fluorophores)
