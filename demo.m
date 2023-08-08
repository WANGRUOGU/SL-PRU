% setup
close all
addpath(genpath(pwd))
warning('off')
filename = 'mixture.tif'; % image file name
image = readMultipageTiff(filename);
load endmembers.mat % import fluorophore spectral signatures

% tuning parameters
% this step can be ignored and directly set the values for lambda1 and
% lambda2
R = size(E, 2); % number of endmembers
ref = cell(1, R);
for r = 1:R % load reference images for tuning
    ref{r} = readMultipageTiff([Fluorophores{r}, '.tif']);
end
lam_seq = @(lam_ub, lam_lb, lam_range) ...
    [0,exp(log(lam_lb):((log(lam_ub)-log(lam_lb))/lam_range):log(lam_ub))];
lam = lam_seq(5, .1, 3);
sz_sam = 1e0;
[l1_PR, l2_PR, prop] = Tuning(ref, E, Fluorophores, lam, lam, sz_sam);

% SL-PRU unmixing
A = UnmixSpIm(image, E, 'SLPRU', l1_PR, l2_PR);

% save the abundances
imwrite(A(:, :, 1), ['./result/Unmixed_', filename]);
for r = 2:R
    imwrite(A(:, :, r), ['./result/Unmixed_', filename], 'WriteMode' , 'append');
end