function [] = peakAlign_prfPCA(thisBA3a_cat,thisBA3b_cat,thisBA1_cat,thisBA2_cat,whichMode,currentModel,noComponents,orthogMe)
%peakAlign_prfPCA

% what if we align to peaks, not just shifting by digit
% 18 is the middle of 35
% See Also prfsummary_baseball_pca_grand.m

if strcmpi(whichMode,'BAA')
    doBAA = 1;
    doPD = 0;
elseif strcmpi(whichMode,'BASETIPS')
    doPD = 1;
    doBAA = 0;
end

%magicBall = 18;
clear x0a x0b x01 x02

for ii = 1:4

    [x0a(ii), y0, v] = find2dMax(thisBA3a_cat(:,:,(ii)));
    [x0b(ii), y0, v] = find2dMax(thisBA3b_cat(:,:,(ii)));
    [x01(ii), y0, v] = find2dMax(thisBA1_cat(:,:,(ii)));
    [x02(ii), y0, v] = find2dMax(thisBA2_cat(:,:,(ii)));

end

GridY = 60; % 20 60
GridX = 105; % 35 105
magicBall = ceil(GridX./2);
GridYM1 = GridY-1;

clear bigman*
bigmana = NaN(GridY,GridX,4);
bigmanb = NaN(GridY,GridX,4);
bigman1 = NaN(GridY,GridX,4);
bigman2 = NaN(GridY,GridX,4);

clear bloop*

% lots of fudging here to make sure that we fill a 20x20 matrix correctly
for ii = 1:4
    bloopa_start(ii) = magicBall-x0a(ii);
    bloopa_end(ii) = magicBall-x0a(ii)+GridYM1;
    if bloopa_start(ii)<0; bloopa_start(ii) = 1; end
    if bloopa_end(ii)>GridX
        bloopa_end(ii) = GridX;
        if length(bloopa_start(ii):bloopa_end(ii)) < GridY
            bloopa_start(ii) = bloopa_end(ii)-GridYM1;
        end
    end
    if length(bloopa_start(ii):bloopa_end(ii))<GridY
        bloopa_end(ii) = bloopa_start(ii)+GridYM1; %very dodgy
    end

    bloopb_start(ii) = magicBall-x0b(ii);
    bloopb_end(ii) = magicBall-x0b(ii)+GridYM1;
    if bloopb_start(ii)<0; bloopb_start(ii) = 1; end
    if bloopb_end(ii)>GridX
        bloopb_end(ii) = GridX;
        if length(bloopb_start(ii):bloopb_end(ii)) < GridY
            bloopb_start(ii) = bloopb_end(ii)-GridYM1;
        end
    end
    if length(bloopb_start(ii):bloopb_end(ii))<GridY
        bloopb_end(ii) = bloopb_start(ii)+GridYM1; %very dodgy
    end


    bloop1_start(ii) = magicBall-x01(ii);
    bloop1_end(ii) = magicBall-x01(ii)+GridYM1;
    if bloop1_start(ii)<0; bloop1_start(ii) = 1; end
    if bloop1_end(ii)>GridX
        bloop1_end(ii) = GridX;
        if length(bloop1_start(ii):bloop1_end(ii)) < GridY
            bloop1_start(ii) = bloop1_end(ii)-GridYM1;
        end
    end
    if length(bloop1_start(ii):bloop1_end(ii))<GridY
        bloop1_end(ii) = bloop1_start(ii)+GridYM1; %very dodgy
    end

    bloop2_start(ii) = magicBall-x02(ii);
    bloop2_end(ii) = magicBall-x02(ii)+GridYM1;
    if bloop2_start(ii)<0; bloop2_start(ii) = 1; end
    if bloop2_end(ii)>GridX
        bloop2_end(ii) = GridX;
        if length(bloop2_start(ii):bloop2_end(ii)) < GridY
            bloop2_start(ii) = bloop2_end(ii)-GridYM1;
        end
    end
    if length(bloop2_start(ii):bloop2_end(ii))<GridY
        bloop2_end(ii) = bloop2_start(ii)+GridYM1; %very dodgy
    end



end


% now fill
for ii = 1:4

    bigmana(:,bloopa_start(ii):bloopa_end(ii),ii) = thisBA3a_cat(:,:,ii);
    bigmanb(:,bloopb_start(ii):bloopb_end(ii),ii) = thisBA3b_cat(:,:,ii);
    bigman1(:,bloop1_start(ii):bloop1_end(ii),ii) = thisBA1_cat(:,:,ii);
    bigman2(:,bloop2_start(ii):bloop2_end(ii),ii) = thisBA2_cat(:,:,ii);

end

bigmean3a_peak  = nanmean(bigmana,3);
bigmean3b_peak  = nanmean(bigmanb,3);
bigmean1_peak  = nanmean(bigman1,3);
bigmean2_peak  = nanmean(bigman2,3);

%%

if doPD

    figure
    tiledlayout(1,4)
    axis square
    colormap jet
    ax1 = nexttile;
    imagesc(ax1,bigmean3a_peak)
    title('PD1')
    ax2 = nexttile;
    imagesc(ax2,bigmean3b_peak)
    title('PD2')
    ax3 = nexttile;
    imagesc(ax3,bigmean1_peak)
    title('PD3')
    ax4 = nexttile;
    imagesc(ax4,bigmean2_peak)
    title('PD4')
    %print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BAA_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')
    if orthogMe ~= 1
        print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BASETIPS_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')
    elseif orthogMe == 1
        print(gcf, sprintf('ORTHOGPEAKalignedaverage_%s_PC%s_BASETIPS_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')
    end
elseif doBAA

    figure
    tiledlayout(1,4)
    axis square
    colormap jet
    ax1 = nexttile;
    imagesc(ax1,bigmean3a_peak)
    title('BA3a')
    ax2 = nexttile;
    imagesc(ax2,bigmean3b_peak)
    title('BA3b')
    ax3 = nexttile;
    imagesc(ax3,bigmean1_peak)
    title('BA1')
    ax4 = nexttile;
    imagesc(ax4,bigmean2_peak)
    title('BA2')
    if orthogMe ~= 1
        print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BAA_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')
    elseif orthogMe == 1
        print(gcf, sprintf('ORTHOGPEAKalignedaverage_%s_PC%s_BAA_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')

    end
    %print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BASETIPS_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')

end






end