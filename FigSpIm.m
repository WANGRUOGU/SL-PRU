function FigSpIm(T, A_NLS, A_SNLS, A_SLNLS, A_SLPRU, Fluorophores)

C = size(T, 3);
% Spectral images by channel
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(3, 8,"TileSpacing","compact")
for c = 1:C
    nexttile
    imshow(T(:, :, c))
    title(['Channel ', num2str(c)],'FontSize',20)
end
% Estimated abundances by method (first half)
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(4, 7,"TileSpacing","compact")
for r = 1:7
    nexttile
    imshow(A_NLS(:, :, r))
    if r == 1
        ylabel('NLS','FontSize',20)
    end
    title(Fluorophores{r},'FontSize',20)
end
for r = 1:7
    nexttile
    imshow(A_SNLS(:, :, r))
    if r == 1
        ylabel('S-NLS','FontSize',20)
    end
end
for r = 1:7
    nexttile
    imshow(A_SLNLS(:, :, r))
    if r == 1
        ylabel('SL-NLS','FontSize',20)
    end
end
for r = 1:7
    nexttile
    imshow(A_SLPRU(:, :, r))
    if r == 1
        ylabel('SL-PRU','FontSize',20)
    end
end
% Estimated abundances by method (second half)
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(4, 6,"TileSpacing","compact")
for r = 8:13
    nexttile
    imshow(A_NLS(:, :, r))
    if r == 8
        ylabel('NLS','FontSize',20)
    end
    title(Fluorophores{r},'FontSize',20)
end
for r = 8:13
    nexttile
    imshow(A_SNLS(:, :, r))
    if r == 8
        ylabel('S-NLS','FontSize',20)
    end
end
for r = 8:13
    nexttile
    imshow(A_SLNLS(:, :, r))
    if r == 8
        ylabel('SL-NLS','FontSize',20)
    end
end
for r = 8:13
    nexttile
    imshow(A_SLPRU(:, :, r))
    if r == 8
        ylabel('SL-PRU','FontSize',20)
    end
end