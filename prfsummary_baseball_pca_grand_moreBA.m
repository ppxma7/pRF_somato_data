%% PCA
% this does the mean of the first PC for each subject, weighted by the
% first PC explained variance
mypath = '/Volumes/nemosine/prfsomato_fits/60x60_moreBA/';
cd(mypath)
close all
clc

%a = 1; % principal component to display
noComponents = 1:3; % how many components to use in grand PCA
% any less than 16 (because we have a 4x4 grid) and it breaks the PCA
thisModel = 3;
doBAA = 1;
doPD = 0;
doMeans = 0;

if doBAA
    if thisModel == 1
        theGRANDpca = load('prf_2ddouble_dd_nov_bigBAA_myPCA_big.mat');
        currentModel = '2d';
    elseif thisModel == 2
        theGRANDpca = load('prf_1d_nodiag2_nov_bigBAA_myPCA_big.mat');
        currentModel = '1d';
    elseif thisModel == 3
        theGRANDpca = load('prf_1dt_nodiag2_nov_bigBAA_myPCA_big.mat');
        currentModel = '1dt';
    elseif thisModel == 4
        theGRANDpca = load('prf_sixteen_nodiag2_com_nov_bigBAA_myPCA_big.mat');
        currentModel = '16';
    end
elseif doPD
    if thisModel == 1
        theGRANDpca = load('prf_2ddouble_dd_nov_bigBASETIPS_myPCA.mat');
        currentModel = '2d';
    elseif thisModel == 2
        theGRANDpca = load('prf_1d_nodiag2_nov_bigBASETIPS_myPCA.mat');
        currentModel = '1d';
    elseif thisModel == 3
        theGRANDpca = load('prf_1dt_nodiag2_nov_bigBASETIPS_myPCA.mat');
        currentModel = '1dt';
    elseif thisModel == 4
        theGRANDpca = load('prf_sixteen_nodiag2_com_nov_bigBASETIPS_myPCA.mat');
        currentModel = '16';
    end
end

if doBAA
    balabel = {'BA3a','BA3b','BA1','BA2','BA4a','BA4p','BA6'};
    splitting = 'BAA';
elseif doPD
    balabel = {'PD1','PD2','PD3','PD4'};
    splitting = 'BASETIPS';
end



%
% because PCA matrices may be different sizes (always 16 rows, because 16
% vars, but can be less than 16 cols, because we had not always 48 voxles
% when we ran first PCA), so need to unpack not into matrix, but into
% another struct

% we want this to be each column is Digit, each row is PD(top is PD2(BA3a,
% bottom is PD5 (BA2

gridX = 60; %20 %4
gridY = gridX;
%gridY = 60; %20 %4
%
superGrid = 105; %105 %35; % this is to do digit align, we need extra space. For the 20x20 case, where each column is 5
% to overlap correctly we need a matrix that has 7 columns,
% 7*5 = 35. If we interp more, so say 60x60, then we still
% need 7 columns, but each col is now 15, so 7*15 = 105.
if superGrid == 35
    chopA = 16:35;
    chopB = 11:30;
    chopC = 6:25;
    chopD = 1:20;
elseif superGrid == 105
    chopA = 46:105;
    chopB = 31:90;
    chopC = 16:75;
    chopD = 1:60;
end


if gridX == 20
    cat1 = 1:20;
    cat2 = 21:40;
    cat3 = 41:60;
    cat4 = 61:80;
elseif gridX == 4
    cat1 = 1:4;
    cat2 = 5:8;
    cat3 = 9:12;
    cat4 = 13:16;
elseif gridX == 60
    cat1 = 1:60;
    cat2 = 61:120;
    cat3 = 121:180;
    cat4 = 181:240;
end




orthogMe = 0; % Which direction do you want to sum, across PD/BAs, then set to 1, 0 = across digs.


%% BEGIN
%
%
    for iSub = 1:8
        for iBA = 1:7
            for iDig = 1:4

                %size error check
                if size(theGRANDpca.myGRANDpca.group{iSub}{iBA,iDig},2) < length(noComponents)
                    noComponents = 1;
                else
                    %thisComponent = noComponents;
                end


                digits.location{iBA,iDig,iSub} = theGRANDpca.myGRANDpca.group{iSub}{iBA,iDig}(:,noComponents);
                digits.explained{iBA,iDig,iSub} = theGRANDpca.myGRANDpca.explained{iSub}{iBA,iDig}(noComponents,:);
            end
        end
    end

    if numel(noComponents) == 1
        noComponentidx = 1;
    else
        noComponentidx = noComponents;
    end

    % weight the first PC by the subjectwise explained variance
    for iSub = 1:8
        for iBA = 1:7
            for iDig = 1:4
                % for one component this is easy
                % for multiple components, we can do a matrix multiplication
                % dot product bla * bla. This is equivalent to
                % sum(bla).*transpose(bla)
                digits.weighted{iBA,iDig,iSub} = digits.location{iBA,iDig,iSub}(:,noComponentidx) * (digits.explained{iBA,iDig,iSub}(noComponentidx)./100);

                % and now put into one stack, each col is a sub
                digits.weightedSubs{iBA,iDig}(:,iSub) = digits.weighted{iBA,iDig,iSub};
                digits.unweightedSubs{iBA,iDig}(:,iSub) = sum(digits.location{iBA,iDig,iSub},2);
            end
        end
    end

    % now find mean weightings per location
    for iBA = 1:7
        for iDig = 1:4
            digits.meanPCAweightedRF{iBA,iDig} = mean(digits.weightedSubs{iBA,iDig},2);
            digits.meanPCAunweightedRF{iBA,iDig} = mean(digits.unweightedSubs{iBA,iDig},2);
            % reshape into RF grid
            digits.meanPCAweightedRFreshaped{iBA,iDig} = reshape(digits.meanPCAweightedRF{iBA,iDig},gridX,gridY);
            digits.meanPCAunweightedRFreshaped{iBA,iDig} = reshape(digits.meanPCAunweightedRF{iBA,iDig},gridX,gridX);
        end
    end

    theDigs = 2:5;

    %% now display
    figure('Position',[100 100 1200 1000])
    tiledlayout(7,4)
    
    for iBa = 1:7
        myBa = balabel{iBa};

        for iDig = 1:4
            myDig = theDigs(iDig);
            nexttile
            imagesc(digits.meanPCAweightedRFreshaped{iBa,iDig});
            [x0, y0, v] = find2dMax(digits.meanPCAweightedRFreshaped{iBa,iDig});
            hold on
            s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
                'markeredgecolor', 'k');
            t_ = text(-2,-8,sprintf('dig:%d, %s, model:%s', myDig, ...
                myBa, currentModel));
            set(t_, 'color', [0 0 0 ], ...
                'fontsize', 12, 'fontname', 'helvetica', 'interpreter', 'none', ...
                'verticalalignment', 'top');
            axis square
            %xticklabels({'D2','D3','D4','D5'})
            %yticklabels({'D2','D3','D4','D5'})
            ax = gca;
            ax.FontSize = 12;
            %caxis([0 0.25])
            colorbar
            colormap jet
        end
        %colorbar

    end
    if doBAA
        print(gcf, sprintf('grandmeanofPCA_weighted_%s_PC%s_BAA_BIG',currentModel, num2str(noComponents)), '-r300', '-dpng')
    elseif doPD
        print(gcf, sprintf('grandmeanofPCA_weighted_%s_PC%s_BASETIPS_BIG',currentModel, num2str(noComponents)), '-r300', '-dpng')
    end


    %% try summing
thisgrid = digits.meanPCAweightedRFreshaped;
if orthogMe ~= 1
    thisBA3a = cell2mat(thisgrid(1,:)); %D2
    thisBA3b = cell2mat(thisgrid(2,:)); %
    thisBA1 = cell2mat(thisgrid(3,:));
    thisBA2 = cell2mat(thisgrid(4,:));
    thisBA4a = cell2mat(thisgrid(5,:));
    thisBA4p = cell2mat(thisgrid(6,:));
    thisBA6 = cell2mat(thisgrid(7,:));
    thisBA3a_cat = cat(3,thisBA3a(:,cat1),thisBA3a(:,cat2),thisBA3a(:,cat3),thisBA3a(:,cat4));
    thisBA3b_cat = cat(3,thisBA3b(:,cat1),thisBA3b(:,cat2),thisBA3b(:,cat3),thisBA3b(:,cat4));
    thisBA1_cat = cat(3,thisBA1(:,cat1),thisBA1(:,cat2),thisBA1(:,cat3),thisBA1(:,cat4));
    thisBA2_cat = cat(3,thisBA2(:,cat1),thisBA2(:,cat2),thisBA2(:,cat3),thisBA2(:,cat4));
    thisBA4a_cat = cat(3,thisBA4a(:,cat1),thisBA4a(:,cat2),thisBA4a(:,cat3),thisBA4a(:,cat4));
    thisBA4p_cat = cat(3,thisBA4p(:,cat1),thisBA4p(:,cat2),thisBA4p(:,cat3),thisBA4p(:,cat4));
    thisBA6_cat = cat(3,thisBA6(:,cat1),thisBA6(:,cat2),thisBA6(:,cat3),thisBA6(:,cat4));

else
    thisBA3a = cell2mat(thisgrid(:,1)); %D2
    thisBA3b = cell2mat(thisgrid(:,2)); %
    thisBA1 = cell2mat(thisgrid(:,3));
    thisBA2 = cell2mat(thisgrid(:,4));
    thisBA4a = cell2mat(thisgrid(:,5));
    thisBA4p = cell2mat(thisgrid(:,6));
    thisBA6 = cell2mat(thisgrid(:,7));
    thisBA3a_cat = cat(3,thisBA3a(cat1,:),thisBA3a(cat2,:),thisBA3a(cat3,:),thisBA3a(cat4,:));
    thisBA3b_cat = cat(3,thisBA3b(cat1,:),thisBA3b(cat2,:),thisBA3b(cat3,:),thisBA3b(cat4,:));
    thisBA1_cat = cat(3,thisBA1(cat1,:),thisBA1(cat2,:),thisBA1(cat3,:),thisBA1(cat4,:));
    thisBA2_cat = cat(3,thisBA2(cat1,:),thisBA2(cat2,:),thisBA2(cat3,:),thisBA2(cat4,:));
    thisBA4a_cat = cat(3,thisBA4a(cat1,:),thisBA4a(cat2,:),thisBA4a(cat3,:),thisBA4a(cat4,:));
    thisBA4p_cat = cat(3,thisBA4p(cat1,:),thisBA4p(cat2,:),thisBA4p(cat3,:),thisBA4p(cat4,:));
    thisBA6_cat = cat(3,thisBA6(cat1,:),thisBA6(cat2,:),thisBA6(cat3,:),thisBA6(cat4,:));

end

thisBA3a_cat_mean = mean(thisBA3a_cat,3);
thisBA3b_cat_mean = mean(thisBA3b_cat,3);
thisBA1_cat_mean = mean(thisBA1_cat,3);
thisBA2_cat_mean = mean(thisBA2_cat,3);
thisBA4a_cat_mean = mean(thisBA4a_cat,3);
thisBA4p_cat_mean = mean(thisBA4p_cat,3);
thisBA6_cat_mean = mean(thisBA6_cat,3);

% plotting
%balabel = {'BA3a','BA3b','BA1','BA2'};
diglabel = {'D2','D3','D4','D5'};
cellcat = {thisBA3a_cat_mean,thisBA3b_cat_mean,...
    thisBA1_cat_mean,thisBA2_cat_mean,...
    thisBA4a_cat_mean,thisBA4p_cat_mean,thisBA6_cat_mean};

figure('Position',[100 100 1600 400])
tiledlayout(2,4)
for figloop = 1:7
    nexttile
    imagesc(cellcat{figloop})
    [x0, y0, v] = find2dMax(cellcat{figloop});
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if orthogMe ~= 1; t_ = text(10,-8,balabel{figloop});else t_ = text(10,-5,diglabel{figloop}); end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet
end

if orthogMe ~= 1
    print(gcf, sprintf('grandmeanofPCA_summed_weighted_%s_PC%s_%s_interped_BIG',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
else
    print(gcf, sprintf('ORTHOGgrandmeanofPCA_summed_weighted_%s_PC%s_%s_interped_BIG',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
end

%% run peak alignment

peakAlign_prfPCA_moreBA(thisBA3a_cat,thisBA3b_cat,thisBA1_cat,thisBA2_cat,thisBA4a_cat, thisBA4p_cat, thisBA6_cat, splitting,currentModel,noComponents)

%% can we rejig things

NANman = NaN(gridX,superGrid,4);

bigman = NANman;
bigman(:,chopA,1) = thisBA3a_cat(:,:,1);
bigman(:,chopB,2) = thisBA3a_cat(:,:,2);
bigman(:,chopC,3) = thisBA3a_cat(:,:,3);
bigman(:,chopD,4) = thisBA3a_cat(:,:,4);
bigmean3a = nanmean(bigman,3);
%figure, imagesc(bigmean3a)

bigman = NANman;
bigman(:,chopA,1) = thisBA3b_cat(:,:,1);
bigman(:,chopB,2) = thisBA3b_cat(:,:,2);
bigman(:,chopC,3) = thisBA3b_cat(:,:,3);
bigman(:,chopD,4) = thisBA3b_cat(:,:,4);
bigmean3b = nanmean(bigman,3);
%figure, imagesc(bigmean3b)

bigman = NANman;
bigman(:,chopA,1) = thisBA1_cat(:,:,1);
bigman(:,chopB,2) = thisBA1_cat(:,:,2);
bigman(:,chopC,3) = thisBA1_cat(:,:,3);
bigman(:,chopD,4) = thisBA1_cat(:,:,4);
bigmean1 = nanmean(bigman,3);
%figure, imagesc(bigmean1)

bigman = NANman;
bigman(:,chopA,1) = thisBA2_cat(:,:,1);
bigman(:,chopB,2) = thisBA2_cat(:,:,2);
bigman(:,chopC,3) = thisBA2_cat(:,:,3);
bigman(:,chopD,4) = thisBA2_cat(:,:,4);
bigmean2 = nanmean(bigman,3);

bigman = NANman;
bigman(:,chopA,1) = thisBA4a_cat(:,:,1);
bigman(:,chopB,2) = thisBA4a_cat(:,:,2);
bigman(:,chopC,3) = thisBA4a_cat(:,:,3);
bigman(:,chopD,4) = thisBA4a_cat(:,:,4);
bigmean4a = nanmean(bigman,3);

bigman = NANman;
bigman(:,chopA,1) = thisBA4p_cat(:,:,1);
bigman(:,chopB,2) = thisBA4p_cat(:,:,2);
bigman(:,chopC,3) = thisBA4p_cat(:,:,3);
bigman(:,chopD,4) = thisBA4p_cat(:,:,4);
bigmean4p = nanmean(bigman,3);

bigman = NANman;
bigman(:,chopA,1) = thisBA6_cat(:,:,1);
bigman(:,chopB,2) = thisBA6_cat(:,:,2);
bigman(:,chopC,3) = thisBA6_cat(:,:,3);
bigman(:,chopD,4) = thisBA6_cat(:,:,4);
bigmean6 = nanmean(bigman,3);

figure
tiledlayout(1,4)
axis square
colormap jet
ax1 = nexttile;
imagesc(ax1,bigmean3a)
title(balabel{1})
ax2 = nexttile;
imagesc(ax2,bigmean3b)
title(balabel{2})
ax3 = nexttile;
imagesc(ax3,bigmean1)
title(balabel{3})
ax4 = nexttile;
imagesc(ax4,bigmean2)
title(balabel{4})
ax5 = nexttile;
imagesc(ax5,bigmean4a)
title(balabel{5})
ax6 = nexttile;
imagesc(ax6,bigmean4p)
title(balabel{6})
ax7 = nexttile;
imagesc(ax7,bigmean6)
title(balabel{7})
if doBAA
    print(gcf, sprintf('alignedaverage_%s_PC%s_%s_interped_BIG',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
elseif doPD
    print(gcf, sprintf('alignedaverage_%s_PC%s_%s_interped_BIG',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
end

%% do means
if doMeans
    mypath = '/Volumes/nemosine/prfsomato_fits';
    cd(mypath)
    close all
    clc

    gridX = 20; %4
    gridY = 20; %4

    doBasetips = 0;
    doBa = 1;



    currentModel = '16';

    if doBa == 1
        sep = 'BAA';
        howmanyrows = 7;
    elseif doBasetips == 1
        sep = 'BASETIPS';
        howmanyrows = 4;
    end

    if strcmpi(currentModel,'16')
        if doBasetips == 1
            theGRANDmean = load('prf_sixteen_nodiag2_com_nov_bigBASETIPS_myMEAN.mat');
        elseif doBa == 1
            theGRANDmean = load('prf_sixteen_nodiag2_com_nov_bigBAA_myMEAN.mat');
        end
    elseif strcmpi(currentModel,'2d')
        if doBasetips == 1
            theGRANDmean = load('prf_2ddouble_dd_novBASETIPS_myMEAN.mat');
        elseif doBa == 1
            theGRANDmean = load('prf_2ddouble_dd_novBAA_myMEAN.mat');
        end
    elseif strcmpi(currentModel,'1d')
        if doBasetips == 1
            theGRANDmean = load('prf_1d_nodiag2_novBASETIPS_myMEAN.mat');
        elseif doBa == 1
            theGRANDmean = load('prf_1d_nodiag2_novBAA_myMEAN.mat');
        end
    elseif strcmpi(currentModel,'1dt')
        if doBasetips == 1
            theGRANDmean = load('prf_1dt_nodiag2_novBASETIPS_myMEAN.mat');
        elseif doBa == 1
            theGRANDmean = load('prf_1dt_nodiag2_novBAA_myMEAN.mat');
        end
    end



    for iSub = 1:8
        for iBA = 1:howmanyrows
            for iDig = 1:4
                digits.location{iBA,iDig,iSub} = theGRANDmean.myGRANDmean.group{iSub}{iBA,iDig};
            end
        end
    end

    for iBA = 1:howmanyrows
        for iDig = 1:4
            digits.means{iBA,iDig} = mean(cell2mat(digits.location(iBA,iDig,:)),3);
        end
    end

    figure('Position',[100 100 1200 1000])
    tiledlayout(howmanyrows,4)
    if doBasetips == 1
        theLocs = {'PD1','PD2','PD3','PD4'};
    elseif doBa == 1
        theLocs = {'BA3a','BA3b','BA1','BA2','BA4a','BA4p','BA6'};
    end
    theDigs = 2:5;
    for iLoc = 1:howmanyrows
        myLoc = theLocs{iLoc};

        for iDig = 1:4
            myDig = theDigs(iDig);
            nexttile
            imagesc(digits.means{iLoc,iDig});
            [x0, y0, v] = find2dMax(digits.means{iLoc,iDig});
            hold on
            s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
                'markeredgecolor', 'k');
            t_ = text(0,-3,sprintf('dig:%d, %s, model:%s', myDig, ...
                myLoc, currentModel));
            set(t_, 'color', [0 0 0 ], ...
                'fontsize', 12, 'fontname', 'helvetica', 'interpreter', 'none', ...
                'verticalalignment', 'top');
            axis square
            ax = gca;
            ax.FontSize = 12;
            colorbar
            colormap jet
        end

    end

    print(gcf, sprintf('gmean_%s_%s',currentModel,sep), '-r300', '-dpng')


    %% try summing
    thisgrid = digits.means;
    thisBA3a = cell2mat(thisgrid(1,:)); %D2
    % thisBA3a_cat = cat(3,thisBA3a(:,1:4),thisBA3a(:,5:8),thisBA3a(:,9:12),thisBA3a(:,13:16));
    % thisBA3a_cat_mean = mean(thisBA3a_cat,3);

    if gridX == 20
        thisBA3a_cat = cat(3,thisBA3a(:,1:20),thisBA3a(:,21:40),thisBA3a(:,41:60),thisBA3a(:,61:80));
    elseif gridX == 4
        thisBA3a_cat = cat(3,thisBA3a(:,1:4),thisBA3a(:,5:8),thisBA3a(:,9:12),thisBA3a(:,13:16));
    end
    thisBA3a_cat_mean = mean(thisBA3a_cat,3);



    % now smooth
    %[thisBA3a_cat_mean_smoothed] = smoothMyRF(thisBA3a_cat_mean);

    thisBA3b = cell2mat(thisgrid(2,:));%D3
    if gridX == 20
        thisBA3b_cat = cat(3,thisBA3b(:,1:20),thisBA3b(:,21:40),thisBA3b(:,41:60),thisBA3b(:,61:80));
    elseif gridX == 4
        thisBA3b_cat = cat(3,thisBA3b(:,1:4),thisBA3b(:,5:8),thisBA3b(:,9:12),thisBA3b(:,13:16));
    end
    thisBA3b_cat_mean = mean(thisBA3b_cat,3);

    %
    % thisBA3b = cell2mat(thisgrid(2,:));%D3
    % thisBA3b_cat = cat(3,thisBA3b(:,1:4),thisBA3b(:,5:8),thisBA3b(:,9:12),thisBA3b(:,13:16));
    % thisBA3b_cat_mean = mean(thisBA3b_cat,3);
    % now smooth
    %[thisBA3b_cat_mean_smoothed] = smoothMyRF(thisBA3b_cat_mean);
    %
    % thisBA1 = cell2mat(thisgrid(3,:));%D4
    % thisBA1_cat = cat(3,thisBA1(:,1:4),thisBA1(:,5:8),thisBA1(:,9:12),thisBA1(:,13:16));
    % thisBA1_cat_mean = mean(thisBA1_cat,3);
    % now smooth
    %[thisBA1_cat_mean_smoothed] = smoothMyRF(thisBA1_cat_mean);
    thisBA1 = cell2mat(thisgrid(3,:));%D4
    if gridX == 20
        thisBA1_cat = cat(3,thisBA1(:,1:20),thisBA1(:,21:40),thisBA1(:,41:60),thisBA1(:,61:80));
    elseif gridX == 4
        thisBA1_cat = cat(3,thisBA1(:,1:4),thisBA1(:,5:8),thisBA1(:,9:12),thisBA1(:,13:16));
    end
    thisBA1_cat_mean = mean(thisBA1_cat,3);

    thisBA2 = cell2mat(thisgrid(4,:));%D4
    if gridX == 20
        thisBA2_cat = cat(3,thisBA2(:,1:20),thisBA2(:,21:40),thisBA2(:,41:60),thisBA2(:,61:80));
    elseif gridX == 4
        thisBA2_cat = cat(3,thisBA2(:,1:4),thisBA2(:,5:8),thisBA2(:,9:12),thisBA2(:,13:16));
    end
    thisBA2_cat_mean = mean(thisBA2_cat,3);

    if doBa == 1

        thisBA4a = cell2mat(thisgrid(5,:));%D4
        if gridX == 20
            thisBA4a_cat = cat(3,thisBA4a(:,1:20),thisBA4a(:,21:40),thisBA4a(:,41:60),thisBA4a(:,61:80));
        elseif gridX == 4
            thisBA4a_cat = cat(3,thisBA4a(:,1:4),thisBA4a(:,5:8),thisBA4a(:,9:12),thisBA4a(:,13:16));
        end
        thisBA4a_cat_mean = mean(thisBA4a_cat,3);

        thisBA4p = cell2mat(thisgrid(6,:));%D4
        if gridX == 20
            thisBA4p_cat = cat(3,thisBA4p(:,1:20),thisBA4p(:,21:40),thisBA4p(:,41:60),thisBA4p(:,61:80));
        elseif gridX == 4
            thisBA4p_cat = cat(3,thisBA4p(:,1:4),thisBA4p(:,5:8),thisBA4p(:,9:12),thisBA4p(:,13:16));
        end
        thisBA4p_cat_mean = mean(thisBA4p_cat,3);

        thisBA6 = cell2mat(thisgrid(7,:));%D4
        if gridX == 20
            thisBA6_cat = cat(3,thisBA6(:,1:20),thisBA6(:,21:40),thisBA6(:,41:60),thisBA6(:,61:80));
        elseif gridX == 4
            thisBA6_cat = cat(3,thisBA6(:,1:4),thisBA6(:,5:8),thisBA6(:,9:12),thisBA6(:,13:16));
        end
        thisBA6_cat_mean = mean(thisBA6_cat,3);

        % thisBA2 = cell2mat(thisgrid(4,:));%D5
        % thisBA2_cat = cat(3,thisBA2(:,1:4),thisBA2(:,5:8),thisBA2(:,9:12),thisBA2(:,13:16));
        % thisBA2_cat_mean = mean(thisBA2_cat,3);
        % now smooth
        %[thisBA2_cat_mean_smoothed] = smoothMyRF(thisBA2_cat_mean);

    end

    figure('Position',[100 100 1600 400])
    if doBa == 1
        tiledlayout(2,4)
    elseif doBasetips == 1
        tiledlayout(1,4)

    end
    nexttile
    imagesc(thisBA3a_cat_mean)
    [x0, y0, v] = find2dMax(thisBA3a_cat_mean);
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD1'); else, t_ = text(10,-4,'BA3a'); end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet

    nexttile
    imagesc(thisBA3b_cat_mean)
    [x0, y0, v] = find2dMax(thisBA3b_cat_mean);
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD2'); else, t_ = text(10,-4,'BA3b'); end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet

    nexttile
    imagesc(thisBA1_cat_mean)
    [x0, y0, v] = find2dMax(thisBA1_cat_mean);
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD3'); else, t_ = text(10,-4,'BA1'); end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet

    nexttile
    imagesc(thisBA2_cat_mean)
    [x0, y0, v] = find2dMax(thisBA2_cat_mean);
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD4'); else, t_ = text(10,-4,'BA2'); end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
        'helvetica', 'interpreter', 'none','verticalalignment','top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet

    if doBa == 1
        nexttile
        imagesc(thisBA4a_cat_mean)
        [x0, y0, v] = find2dMax(thisBA4a_cat_mean);
        hold on
        s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
            'markeredgecolor', 'k');
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD4'); else, t_ = text(10,-4,'BA4a'); end
        set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
            'helvetica', 'interpreter', 'none','verticalalignment','top');
        ax = gca;
        ax.FontSize = 12;
        colorbar
        axis square
        colormap jet

        nexttile
        imagesc(thisBA4p_cat_mean)
        [x0, y0, v] = find2dMax(thisBA4p_cat_mean);
        hold on
        s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
            'markeredgecolor', 'k');
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD4'); else, t_ = text(10,-4,'BA4p'); end
        set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
            'helvetica', 'interpreter', 'none','verticalalignment','top');
        ax = gca;
        ax.FontSize = 12;
        colorbar
        axis square
        colormap jet

        nexttile
        imagesc(thisBA6_cat_mean)
        [x0, y0, v] = find2dMax(thisBA6_cat_mean);
        hold on
        s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
            'markeredgecolor', 'k');
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD4'); else, t_ = text(10,-4,'BA6'); end
        set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
            'helvetica', 'interpreter', 'none','verticalalignment','top');
        ax = gca;
        ax.FontSize = 12;
        colorbar
        axis square
        colormap jet
    end

    if strcmpi(sep,'BASETIPS')
        print(gcf, sprintf('gsummean_%s_%s_BASETIPS_raw',currentModel, sep), '-r300', '-dpng')
    elseif strcmpi(sep,'BAA')
        print(gcf, sprintf('gsummean_%s_%s_BAA_raw_BIG',currentModel, sep), '-r300', '-dpng')
    end
end
% INTERP
% figure('Position',[100 100 1600 400])
% tiledlayout(1,4)
% nexttile
% imagesc(thisBA3a_cat_mean_smoothed)
% [x0, y0, v] = find2dMax(thisBA3a_cat_mean);
% hold on
% s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%     'markeredgecolor', 'k');
% if strcmpi(sep,'BASETIPS'); t_ = text(10,-1,'PD1'); else, t_ = text(10,-1,'BA3a'); end
% set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
% ax = gca;
% ax.FontSize = 12;
% colorbar
% axis square
%
%
% nexttile
% imagesc(thisBA3b_cat_mean_smoothed)
% [x0, y0, v] = find2dMax(thisBA3b_cat_mean);
% hold on
% s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%     'markeredgecolor', 'k');
% if strcmpi(sep,'BASETIPS'); t_ = text(10,-1,'PD2'); else, t_ = text(10,-1,'BA3b'); end
% set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
% ax = gca;
% ax.FontSize = 12;
% colorbar
% axis square
%
%
% nexttile
% imagesc(thisBA1_cat_mean_smoothed)
% [x0, y0, v] = find2dMax(thisBA1_cat_mean);
% hold on
% s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%     'markeredgecolor', 'k');
% if strcmpi(sep,'BASETIPS'); t_ = text(10,-1,'PD3'); else, t_ = text(10,-1,'BA1'); end
% set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
% ax = gca;
% ax.FontSize = 12;
% colorbar
% axis square
%
% nexttile
% imagesc(thisBA2_cat_mean_smoothed)
% [x0, y0, v] = find2dMax(thisBA2_cat_mean);
% hold on
% s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%     'markeredgecolor', 'k');
% if strcmpi(sep,'BASETIPS'); t_ = text(10,-1,'PD4'); else, t_ = text(10,-1,'BA2'); end
% set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
%     'helvetica', 'interpreter', 'none','verticalalignment','top');
% ax = gca;
% ax.FontSize = 12;
% colorbar
% axis square
%
% if strcmpi(sep,'BASETIPS')
%     print(gcf, sprintf('gsummean_%s_%s_BASETIPS_interp',currentModel, sep), '-r300', '-dpng')
% elseif strcmpi(sep,'BAA')
%     print(gcf, sprintf('gsummean_%s_%s_BAA_interp',currentModel, sep), '-r300', '-dpng')
%
% end
%




