%% PCA
% this does the mean of the first PC for each subject, weighted by the
% first PC explained variance
close all
close all hidden
clear all

mypath = '/Volumes/nemosine/prfsomato_fits/60x60/';
cd(mypath)

clc

%a = 1; % principal component to display
noComponents = 1:3; % how many components to use in grand PCA
% any less than 16 (because we have a 4x4 grid) and it breaks the PCA

%model 1 = 2d
%model 2 = 1d
%model 3 = 1dt
%model 4 = 16

thisModel = 3;
doBAA = 0;
doPD = 1;
doMeans = 0;
orthogMe = 1; % Which direction do you want to sum, across PD/BAs, then set to 1, 0 = across digs.

if doBAA
    if thisModel == 1
        theGRANDpca = load('prf_2ddouble_dd_novBAA_myPCA.mat');
        currentModel = '2d';
    elseif thisModel == 2
        theGRANDpca = load('prf_1d_nodiag2_novBAA_myPCA.mat');
        currentModel = '1d';
    elseif thisModel == 3
        theGRANDpca = load('prf_1dt_nodiag2_novBAA_myPCA.mat');
        currentModel = '1dt';
    elseif thisModel == 4
        theGRANDpca = load('prf_sixteen_nodiag2_com_novBAA_myPCA.mat');
        currentModel = '16';
    end
elseif doPD
    if thisModel == 1
        theGRANDpca = load('prf_2ddouble_dd_novBASETIPS_myPCA.mat');
        currentModel = '2d';
    elseif thisModel == 2
        theGRANDpca = load('prf_1d_nodiag2_novBASETIPS_myPCA.mat');
        currentModel = '1d';
    elseif thisModel == 3
        theGRANDpca = load('prf_1dt_nodiag2_novBASETIPS_myPCA.mat');
        currentModel = '1dt';
    elseif thisModel == 4
        theGRANDpca = load('prf_sixteen_nodiag2_com_novBASETIPS_myPCA.mat');
        currentModel = '16';
    end
end


if doBAA
    balabel = {'BA3a','BA3b','BA1','BA2'};
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







%% BEGIN
for iSub = 1:10
    for iBA = 1:4
        for iDig = 1:4


            %size error check
            if size(theGRANDpca.myGRANDpca.group{iSub}{iBA,iDig},2) < length(noComponents)
                noComponents = 1;
            else
                %noComponents = noComponents;
            end
            

%             digits.location{iBA,iDig,iSub} = theGRANDpca.myGRANDpca.group{iSub}{iBA,iDig}(:,thisComponent);
%             digits.explained{iBA,iDig,iSub} = theGRANDpca.myGRANDpca.explained{iSub}{iBA,iDig}(thisComponent,:);


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
for iSub = 1:10
    for iBA = 1:4
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
for iBA = 1:4
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
tiledlayout(4,4)
for iBa = 1:4
    myBa = balabel{iBa};

    for iDig = 1:4
        myDig = theDigs(iDig);
        nexttile
        imagesc(digits.meanPCAweightedRFreshaped{iBa,iDig});
        [x0, y0, v] = find2dMax(digits.meanPCAweightedRFreshaped{iBa,iDig});
        hold on
        s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
            'markeredgecolor', 'k');
        t_ = text(-2,-5,sprintf('dig:%d, %s, model:%s', myDig, ...
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
    print(gcf, sprintf('grandmeanofPCA_weighted_%s_PC%s_BAA',currentModel, num2str(noComponents)), '-r300', '-dpng')
elseif doPD
    print(gcf, sprintf('grandmeanofPCA_weighted_%s_PC%s_BASETIPS',currentModel, num2str(noComponents)), '-r300', '-dpng')
end


%% try summing
thisgrid = digits.meanPCAweightedRFreshaped;
if orthogMe ~= 1
    thisBA3a = cell2mat(thisgrid(1,:)); %D2
    thisBA3b = cell2mat(thisgrid(2,:)); %
    thisBA1 = cell2mat(thisgrid(3,:));
    thisBA2 = cell2mat(thisgrid(4,:));
    thisBA3a_cat = cat(3,thisBA3a(:,cat1),thisBA3a(:,cat2),thisBA3a(:,cat3),thisBA3a(:,cat4));
    thisBA3b_cat = cat(3,thisBA3b(:,cat1),thisBA3b(:,cat2),thisBA3b(:,cat3),thisBA3b(:,cat4));
    thisBA1_cat = cat(3,thisBA1(:,cat1),thisBA1(:,cat2),thisBA1(:,cat3),thisBA1(:,cat4));
    thisBA2_cat = cat(3,thisBA2(:,cat1),thisBA2(:,cat2),thisBA2(:,cat3),thisBA2(:,cat4));
else
    thisBA3a = cell2mat(thisgrid(:,1)); %D2
    thisBA3b = cell2mat(thisgrid(:,2)); %
    thisBA1 = cell2mat(thisgrid(:,3));
    thisBA2 = cell2mat(thisgrid(:,4));
    thisBA3a_cat = cat(3,thisBA3a(cat1,:),thisBA3a(cat2,:),thisBA3a(cat3,:),thisBA3a(cat4,:));
    thisBA3b_cat = cat(3,thisBA3b(cat1,:),thisBA3b(cat2,:),thisBA3b(cat3,:),thisBA3b(cat4,:));
    thisBA1_cat = cat(3,thisBA1(cat1,:),thisBA1(cat2,:),thisBA1(cat3,:),thisBA1(cat4,:));
    thisBA2_cat = cat(3,thisBA2(cat1,:),thisBA2(cat2,:),thisBA2(cat3,:),thisBA2(cat4,:));
end

thisBA3a_cat_mean = mean(thisBA3a_cat,3);
thisBA3b_cat_mean = mean(thisBA3b_cat,3);
thisBA1_cat_mean = mean(thisBA1_cat,3);
thisBA2_cat_mean = mean(thisBA2_cat,3);



% plotting
%balabel = {'BA3a','BA3b','BA1','BA2'};
diglabel = {'D2','D3','D4','D5'};
cellcat = {thisBA3a_cat_mean,thisBA3b_cat_mean,thisBA1_cat_mean,thisBA2_cat_mean};

figure('Position',[100 100 1600 400])
tiledlayout(1,4)
for figloop = 1:4
    nexttile
    imagesc(cellcat{figloop})
    [x0, y0, v] = find2dMax(cellcat{figloop});
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if orthogMe ~= 1; t_ = text(10,-5,balabel{figloop});else t_ = text(10,-5,diglabel{figloop}); end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', 'helvetica', 'interpreter', 'none','verticalalignment', 'top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet
end

if orthogMe ~= 1
    print(gcf, sprintf('grandmeanofPCA_summed_weighted_%s_PC%s_%s_interped',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
else
    print(gcf, sprintf('ORTHOGgrandmeanofPCA_summed_weighted_%s_PC%s_%s_interped',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
end

% run peak alignment

peakAlign_prfPCA(thisBA3a_cat,thisBA3b_cat,thisBA1_cat,thisBA2_cat,splitting,currentModel,noComponents,orthogMe)

% can we rejig things

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
if doBAA
    print(gcf, sprintf('alignedaverage_%s_PC%s_%s_interped',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
elseif doPD
    print(gcf, sprintf('alignedaverage_%s_PC%s_%s_interped',currentModel, num2str(noComponents),splitting), '-r300', '-dpng')
end


%% do means
if doMeans

    orthogMe = 1;
    mypath = '/Volumes/nemosine/prfsomato_fits';
    cd(mypath)
    close all
    clc
    gridX = 20; %4
    gridY = 20; %4

    doBasetips = 1;
    doBa = 0;

    currentModel = '16';

    if doBa == 1
        sep = 'BAA';
    elseif doBasetips == 1
        sep = 'BASETIPS';
    end

    if strcmpi(currentModel,'16')
        if doBasetips == 1
            theGRANDmean = load('prf_sixteen_nodiag2_com_novBASETIPS_myMEAN.mat');
        elseif doBa == 1
            theGRANDmean = load('prf_sixteen_nodiag2_com_novBAA_myMEAN.mat');
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



    for iSub = 1:10
        for iBA = 1:4
            for iDig = 1:4
                digits.location{iBA,iDig,iSub} = theGRANDmean.myGRANDmean.group{iSub}{iBA,iDig};
            end
        end
    end

    for iBA = 1:4
        for iDig = 1:4
            digits.means{iBA,iDig} = mean(cell2mat(digits.location(iBA,iDig,:)),3);
        end
    end

    figure('Position',[100 100 1200 1000])
    tiledlayout(4,4)
    if doBasetips == 1
        theLocs = {'PD1','PD2','PD3','PD4'};
    elseif doBa == 1
        theLocs = {'BA3a','BA3b','BA1','BA2'};
    end
    theDigs = 2:5;
    for iLoc = 1:4
        myLoc = theLocs{iLoc};

        for iDig = 1:4
            myDig = theDigs(iDig);
            nexttile
            imagesc(digits.means{iLoc,iDig});
            [x0, y0, v] = find2dMax(digits.means{iLoc,iDig});
            hold on
            s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
                'markeredgecolor', 'k');
            t_ = text(0,-2,sprintf('dig:%d, %s, model:%s', myDig, ...
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


    % try summing
    thisgrid = digits.means;

    if orthogMe ~= 1; thisBA3a = cell2mat(thisgrid(1,:)); else; thisBA3a = cell2mat(thisgrid(:,1)); end
    if gridX == 20
        if orthogMe ~= 1
            thisBA3a_cat = cat(3,thisBA3a(:,1:20),thisBA3a(:,21:40),thisBA3a(:,41:60),thisBA3a(:,61:80));
        else
            thisBA3a_cat = cat(3,thisBA3a(1:20,:),thisBA3a(21:40,:),thisBA3a(41:60,:),thisBA3a(61:80,:));
        end
    elseif gridX == 4
        thisBA3a_cat = cat(3,thisBA3a(:,1:4),thisBA3a(:,5:8),thisBA3a(:,9:12),thisBA3a(:,13:16));
    end
    thisBA3a_cat_mean = mean(thisBA3a_cat,3);

    if orthogMe ~= 1; thisBA3b = cell2mat(thisgrid(2,:)); else; thisBA3b = cell2mat(thisgrid(:,2)); end
    if gridX == 20
        if orthogMe ~= 1
            thisBA3b_cat = cat(3,thisBA3b(:,1:20),thisBA3b(:,21:40),thisBA3b(:,41:60),thisBA3b(:,61:80));
        else
            thisBA3b_cat = cat(3,thisBA3b(1:20,:),thisBA3b(21:40,:),thisBA3b(41:60,:),thisBA3b(61:80,:));
        end
    elseif gridX == 4
        thisBA3b_cat = cat(3,thisBA3b(:,1:4),thisBA3b(:,5:8),thisBA3b(:,9:12),thisBA3b(:,13:16));
    end
    thisBA3b_cat_mean = mean(thisBA3b_cat,3);

    if orthogMe ~= 1; thisBA1 = cell2mat(thisgrid(3,:)); else; thisBA1 = cell2mat(thisgrid(:,3)); end
    if gridX == 20
        if orthogMe ~= 1
            thisBA1_cat = cat(3,thisBA1(:,1:20),thisBA1(:,21:40),thisBA1(:,41:60),thisBA1(:,61:80));
        else
            thisBA1_cat = cat(3,thisBA1(1:20,:),thisBA1(21:40,:),thisBA1(41:60,:),thisBA1(61:80,:));
        end
    elseif gridX == 4
        thisBA1_cat = cat(3,thisBA1(:,1:4),thisBA1(:,5:8),thisBA1(:,9:12),thisBA1(:,13:16));
    end
    thisBA1_cat_mean = mean(thisBA1_cat,3);

    if orthogMe ~= 1; thisBA2 = cell2mat(thisgrid(4,:)); else; thisBA2 = cell2mat(thisgrid(:,4)); end
    if gridX == 20
        if orthogMe ~= 1
            thisBA2_cat = cat(3,thisBA2(:,1:20),thisBA2(:,21:40),thisBA2(:,41:60),thisBA2(:,61:80));
        else
            thisBA2_cat = cat(3,thisBA2(1:20,:),thisBA2(21:40,:),thisBA2(41:60,:),thisBA2(61:80,:));
        end
    elseif gridX == 4
        thisBA2_cat = cat(3,thisBA2(:,1:4),thisBA2(:,5:8),thisBA2(:,9:12),thisBA2(:,13:16));
    end
    thisBA2_cat_mean = mean(thisBA2_cat,3);


    figure('Position',[100 100 1600 400])
    tiledlayout(1,4)
    nexttile
    imagesc(thisBA3a_cat_mean)
    [x0, y0, v] = find2dMax(thisBA3a_cat_mean);
    hold on
    s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
        'markeredgecolor', 'k');
    if orthogMe == 1; t_ = text(10,-1,'D2'); elseif orthogMe ~= 1
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD1'); elseif strcmpi(sep,'BAA'); t_ = text(10,-1,'BA3a'); end
    end
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
    %if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD2'); else, t_ = text(10,-1,'BA3a'); end
    if orthogMe == 1; t_ = text(10,-1,'D3'); elseif orthogMe ~= 1
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD2'); elseif strcmpi(sep,'BAA'); t_ = text(10,-1,'BA3b'); end
    end
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
    %if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD3'); else, t_ = text(10,-1,'BA3a'); end
    if orthogMe == 1; t_ = text(10,-1,'D4'); elseif orthogMe ~= 1
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD3'); elseif strcmpi(sep,'BAA'); t_ = text(10,-1,'BA1'); end
    end
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
    %if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD4'); else, t_ = text(10,-1,'BA3a'); end
    if orthogMe == 1; t_ = text(10,-1,'D5'); elseif orthogMe ~= 1
        if strcmpi(sep,'BASETIPS'); t_ = text(2,0,'PD4'); elseif strcmpi(sep,'BAA'); t_ = text(10,-1,'BA2'); end
    end
    set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
        'helvetica', 'interpreter', 'none','verticalalignment','top');
    ax = gca;
    ax.FontSize = 12;
    colorbar
    axis square
    colormap jet

    if orthogMe == 1
        print(gcf, sprintf('ORTHOGgsummean_%s_%s_interp',currentModel, sep), '-r300', '-dpng')
    elseif orthogMe ~= 1
        if strcmpi(sep,'BASETIPS')
            print(gcf, sprintf('gsummean_%s_%s_BASETIPS_interp',currentModel, sep), '-r300', '-dpng')
        elseif strcmpi(sep,'BAA')
            print(gcf, sprintf('gsummean_%s_%s_BAA_interp',currentModel, sep), '-r300', '-dpng')
        end
    end


    %% can we rejig things
    bigman = NaN(20,35,4);
    bigman(:,16:35,1) = thisBA3a_cat(:,:,1);
    bigman(:,11:30,2) = thisBA3a_cat(:,:,2);
    bigman(:,6:25,3) = thisBA3a_cat(:,:,3);
    bigman(:,1:20,4) = thisBA3a_cat(:,:,4);
    bigmean3a = nanmean(bigman,3);
    %figure, imagesc(bigmean3a)

    bigman = NaN(20,35,4);
    bigman(:,16:35,1) = thisBA3b_cat(:,:,1);
    bigman(:,11:30,2) = thisBA3b_cat(:,:,2);
    bigman(:,6:25,3) = thisBA3b_cat(:,:,3);
    bigman(:,1:20,4) = thisBA3b_cat(:,:,4);
    bigmean3b = nanmean(bigman,3);
    %figure, imagesc(bigmean3b)

    bigman = NaN(20,35,4);
    bigman(:,16:35,1) = thisBA1_cat(:,:,1);
    bigman(:,11:30,2) = thisBA1_cat(:,:,2);
    bigman(:,6:25,3) = thisBA1_cat(:,:,3);
    bigman(:,1:20,4) = thisBA1_cat(:,:,4);
    bigmean1 = nanmean(bigman,3);
    %figure, imagesc(bigmean1)

    bigman = NaN(20,35,4);
    bigman(:,16:35,1) = thisBA2_cat(:,:,1);
    bigman(:,11:30,2) = thisBA2_cat(:,:,2);
    bigman(:,6:25,3) = thisBA2_cat(:,:,3);
    bigman(:,1:20,4) = thisBA2_cat(:,:,4);
    bigmean2 = nanmean(bigman,3);
    %figure, imagesc(bigmean2)

    if strcmpi(sep,'BASETIPS')
        figure
        tiledlayout(1,4)
        axis square
        colormap jet
        ax1 = nexttile;
        imagesc(ax1,bigmean3a)
        title('PD1')
        ax2 = nexttile;
        imagesc(ax2,bigmean3b)
        title('PD2')
        ax3 = nexttile;
        imagesc(ax3,bigmean1)
        title('PD3')
        ax4 = nexttile;
        imagesc(ax4,bigmean2)
        title('PD4')
        print(gcf, sprintf('alignedaverage_%s_%s_BASETIPS_interped',currentModel, sep), '-r300', '-dpng')
        %print(gcf, sprintf('gsummean_%s_%s_BASETIPS_interp',currentModel, sep), '-r300', '-dpng')
    elseif strcmpi(sep,'BAA')
        figure
        tiledlayout(1,4)
        axis square
        colormap jet
        ax1 = nexttile;
        imagesc(ax1,bigmean3a)
        title('BA3a')
        ax2 = nexttile;
        imagesc(ax2,bigmean3b)
        title('BA3b')
        ax3 = nexttile;
        imagesc(ax3,bigmean1)
        title('BA1')
        ax4 = nexttile;
        imagesc(ax4,bigmean2)
        title('BA2')
        print(gcf, sprintf('alignedaverage_%s_%s_BAA_interped',currentModel, sep), '-r300', '-dpng')
        %print(gcf, sprintf('gsummean_%s_%s_BAA_interp',currentModel, sep), '-r300', '-dpng')
    end





end





