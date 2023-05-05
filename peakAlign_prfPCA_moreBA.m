function [] = peakAlign_prfPCA_moreBA(thisBA3a_cat,thisBA3b_cat,thisBA1_cat,thisBA2_cat,thisBA4a_cat,thisBA4p_cat,thisBA6_cat,whichMode,currentModel,noComponents)
%peakAlign_prfPCA

% what if we align to peaks, not just shifting by digit
% 18 is the middle of 35
% See Also prfsummary_baseball_pca_grand.m

if strcmpi(whichMode,'BAA')
    doBAA = 1;
    doPD = 0;
elseif strcmpi(whichMode,'PD')
    doPD = 1;
    doBAA = 0;
end

%magicBall = 18;
clear x0a x0b x01 x02 x04a x04p x06

for ii = 1:4

    [x0a(ii), y0, v] = find2dMax(thisBA3a_cat(:,:,(ii)));
    [x0b(ii), y0, v] = find2dMax(thisBA3b_cat(:,:,(ii)));
    [x01(ii), y0, v] = find2dMax(thisBA1_cat(:,:,(ii)));
    [x02(ii), y0, v] = find2dMax(thisBA2_cat(:,:,(ii)));
    [x04a(ii), y0, v] = find2dMax(thisBA4a_cat(:,:,(ii)));
    [x04p(ii), y0, v] = find2dMax(thisBA4p_cat(:,:,(ii)));
    [x06(ii), y0, v] = find2dMax(thisBA6_cat(:,:,(ii)));

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
bigman4a = NaN(GridY,GridX,4);
bigman4p = NaN(GridY,GridX,4);
bigman6 = NaN(GridY,GridX,4);

clear bloop*

% lots of fudging here to make sure that we fill a 20x20 matrix correctly
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

    bloop4a_start(ii) = magicBall-x04a(ii);
    bloop4a_end(ii) = magicBall-x04a(ii)+GridYM1;
    if bloop4a_start(ii)<0; bloop4a_start(ii) = 1; end
    if bloop4a_end(ii)>GridX
        bloop4a_end(ii) = GridX;
        if length(bloop4a_start(ii):bloop4a_end(ii)) < GridY
            bloop4a_start(ii) = bloop4a_end(ii)-GridYM1;
        end
    end
    if length(bloop4a_start(ii):bloop4a_end(ii))<GridY
        bloop4a_end(ii) = bloop4a_start(ii)+GridYM1; %very dodgy
    end

    bloop4p_start(ii) = magicBall-x04p(ii);
    bloop4p_end(ii) = magicBall-x04p(ii)+GridYM1;
    if bloop4p_start(ii)<0; bloop4p_start(ii) = 1; end
    if bloop4p_end(ii)>GridX
        bloop4p_end(ii) = GridX;
        if length(bloop4p_start(ii):bloop4p_end(ii)) < GridY
            bloop4p_start(ii) = bloop4p_end(ii)-GridYM1;
        end
    end
    if length(bloop4p_start(ii):bloop4p_end(ii))<GridY
        bloop4p_end(ii) = bloop4p_start(ii)+GridYM1; %very dodgy
    end

    bloop6_start(ii) = magicBall-x06(ii);
    bloop6_end(ii) = magicBall-x06(ii)+GridYM1;
    if bloop6_start(ii)<0; bloop6_start(ii) = 1; end
    if bloop6_end(ii)>GridX
        bloop4a_end(ii) = GridX;
        if length(bloop6_start(ii):bloop6_end(ii)) < GridY
            bloop6_start(ii) = bloop6_end(ii)-GridYM1;
        end
    end
    if length(bloop6_start(ii):bloop6_end(ii))<GridY
        bloop6_end(ii) = bloop6_start(ii)+GridYM1; %very dodgy
    end



end
% 
% bloopa_end-bloopa_start
% bloopb_end-bloopb_start
% bloop1_end-bloop1_start
% bloop2_end-bloop2_start
% bloop4a_end-bloop4a_start
% bloop4p_end-bloop4p_start
% bloop6_end-bloop6_start

% now fill
for ii = 1:4

    bigmana(:,bloopa_start(ii):bloopa_end(ii),ii) = thisBA3a_cat(:,:,ii);
    bigmanb(:,bloopb_start(ii):bloopb_end(ii),ii) = thisBA3b_cat(:,:,ii);
    bigman1(:,bloop1_start(ii):bloop1_end(ii),ii) = thisBA1_cat(:,:,ii);
    bigman2(:,bloop2_start(ii):bloop2_end(ii),ii) = thisBA2_cat(:,:,ii);
    bigman4a(:,bloop4a_start(ii):bloop4a_end(ii),ii) = thisBA4a_cat(:,:,ii);
    bigman4p(:,bloop4p_start(ii):bloop4p_end(ii),ii) = thisBA4p_cat(:,:,ii);
    bigman6(:,bloop6_start(ii):bloop6_end(ii),ii) = thisBA6_cat(:,:,ii);
end

bigmean3a_peak  = nanmean(bigmana,3);
bigmean3b_peak  = nanmean(bigmanb,3);
bigmean1_peak  = nanmean(bigman1,3);
bigmean2_peak  = nanmean(bigman2,3);
bigmean4a_peak  = nanmean(bigman4a,3);
bigmean4p_peak  = nanmean(bigman4p,3);
bigmean6_peak  = nanmean(bigman6,3);

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
%     ax5 = nexttile;
%     imagesc(ax5,bigmean4a_peak)
%     title('PD5')
%     ax6 = nexttile;
%     imagesc(ax6,bigmean4p_peak)
%     title('PD6')
%     ax7 = nexttile;
%     imagesc(ax7,bigmean6_peak)
%     title('PD7')
    %print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BAA_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')
    print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BASETIPS_interped_BIG',currentModel, num2str(noComponents)), '-r300', '-dpng')
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
    ax5 = nexttile;
    imagesc(ax5,bigmean4a_peak)
    title('BA4a')
    ax6 = nexttile;
    imagesc(ax6,bigmean4p_peak)
    title('BA4p')
    ax7 = nexttile;
    imagesc(ax7,bigmean6_peak)
    title('BA6')
    print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BAA_interped_BIG',currentModel, num2str(noComponents)), '-r300', '-dpng')
    %print(gcf, sprintf('PEAKalignedaverage_%s_PC%s_BASETIPS_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')

end






end