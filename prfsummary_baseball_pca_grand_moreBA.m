%% PCA
% this does the mean of the first PC for each subject, weighted by the
% first PC explained variance
mypath='pRF_somato_data/data/';
cd(mypath)
close all
clc

%a = 1; % principal component to display
noComponents = 1:3; % how many components to use in grand PCA
% any less than 16 (because we have a 4x4 grid) and it breaks the PCA

%
%theGRANDpca_basetips = load('prf_1d_nodiag2_nov_bigBASETIPS_myPCA.mat');
%theGRANDpca = load('prf_1d_nodiag2_nov_bigBAA_myPCA.mat');

%theGRANDpca_basetips = load('prf_1dt_nodiag2_nov_bigBASETIPS_myPCA.mat');
%theGRANDpca = load('prf_1dt_nodiag2_nov_bigBAA_myPCA.mat');

theGRANDpca_basetips = load('prf_sixteen_nodiag2_com_nov_bigBASETIPS_myPCA.mat');
theGRANDpca = load('prf_sixteen_nodiag2_com_nov_bigBAA_myPCA.mat');
% 
% theGRANDpca_basetips = load('prf_2ddouble_dd_nov_bigBASETIPS_myPCA.mat');
% theGRANDpca = load('prf_2ddouble_dd_nov_bigBAA_myPCA.mat');

currentModel = '16';
%
% because PCA matrices may be different sizes (always 16 rows, because 16
% vars, but can be less than 16 cols, because we had not always 48 voxles
% when we ran first PCA), so need to unpack not into matrix, but into
% another struct

% we want this to be each column is Digit, each row is PD(top is PD2(BA3a,
% bottom is PD5 (BA2

gridX = 20; %4
gridY = 20; %4
%
%% FIRST BA
%
%
%
for iSub = 1:8
    for iBA = 1:7
        for iDig = 1:4
            
            %size error check
            if size(theGRANDpca.myGRANDpca.group{iSub}{iBA,iDig},2) < length(noComponents)
                thisComponent = 1;
            else
                thisComponent = noComponents;
            end
            
            
            digits.location{iBA,iDig,iSub} = theGRANDpca.myGRANDpca.group{iSub}{iBA,iDig}(:,thisComponent);
            digits.explained{iBA,iDig,iSub} = theGRANDpca.myGRANDpca.explained{iSub}{iBA,iDig}(thisComponent,:);
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
            %size error check
            if size(digits.location{iBA,iDig,iSub},2) < length(noComponentidx)
                thisComponent = 1;
            else
                thisComponent = noComponentidx;
            end
            
            % for one component this is easy
            % for multiple components, we can do a matrix multiplication
            % dot product bla * bla. This is equivalent to
            % sum(bla).*transpose(bla)
            digits.weighted{iBA,iDig,iSub} = digits.location{iBA,iDig,iSub}(:,thisComponent) * (digits.explained{iBA,iDig,iSub}(thisComponent)./100);
            
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
theBaa = {'BA3a','BA3b','BA1','BA2','BA4a','BA4p','BA6'};
for iBa = 1:7
    myBa = theBaa{iBa};
    
    for iDig = 1:4
        myDig = theDigs(iDig);
        nexttile
        imagesc(digits.meanPCAweightedRFreshaped{iBa,iDig});
        [x0, y0, v] = find2dMax(digits.meanPCAweightedRFreshaped{iBa,iDig});
        hold on
        s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
            'markeredgecolor', 'k');
        t_ = text(-2,-3,sprintf('dig:%d, %s, model:%s', myDig, ...
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
print(gcf, sprintf('grandmeanofPCA_weighted_%s_PC%s_BAA_BIG',currentModel, num2str(noComponents)), '-r300', '-dpng')

%% try summing
thisgrid = digits.meanPCAweightedRFreshaped;
thisBA3a = cell2mat(thisgrid(1,:)); %D2
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
% now smooth
%[thisBA3b_cat_mean_smoothed] = smoothMyRF(thisBA3b_cat_mean);

thisBA1 = cell2mat(thisgrid(3,:));%D4
if gridX == 20
    thisBA1_cat = cat(3,thisBA1(:,1:20),thisBA1(:,21:40),thisBA1(:,41:60),thisBA1(:,61:80));
elseif gridX == 4
    thisBA1_cat = cat(3,thisBA1(:,1:4),thisBA1(:,5:8),thisBA1(:,9:12),thisBA1(:,13:16));
end
thisBA1_cat_mean = mean(thisBA1_cat,3);
% now smooth
%[thisBA1_cat_mean_smoothed] = smoothMyRF(thisBA1_cat_mean);

thisBA2 = cell2mat(thisgrid(4,:));%D5
if gridX == 20
    thisBA2_cat = cat(3,thisBA2(:,1:20),thisBA2(:,21:40),thisBA2(:,41:60),thisBA2(:,61:80));
elseif gridX == 4
    thisBA2_cat = cat(3,thisBA2(:,1:4),thisBA2(:,5:8),thisBA2(:,9:12),thisBA2(:,13:16));
end
thisBA2_cat_mean = mean(thisBA2_cat,3);

thisBA4a = cell2mat(thisgrid(5,:));%D5
if gridX == 20
    thisBA4a_cat = cat(3,thisBA4a(:,1:20),thisBA4a(:,21:40),thisBA4a(:,41:60),thisBA4a(:,61:80));
elseif gridX == 4
    thisBA4a_cat = cat(3,thisBA4a(:,1:4),thisBA4a(:,5:8),thisBA4a(:,9:12),thisBA4a(:,13:16));
end
thisBA4a_cat_mean = mean(thisBA4a_cat,3);

thisBA4p = cell2mat(thisgrid(6,:));%D5
if gridX == 20
    thisBA4p_cat = cat(3,thisBA4p(:,1:20),thisBA4p(:,21:40),thisBA4p(:,41:60),thisBA4p(:,61:80));
elseif gridX == 4
    thisBA4p_cat = cat(3,thisBA4p(:,1:4),thisBA4p(:,5:8),thisBA4p(:,9:12),thisBA4p(:,13:16));
end
thisBA4p_cat_mean = mean(thisBA4p_cat,3);

thisBA6 = cell2mat(thisgrid(7,:));%D5
if gridX == 20
    thisBA6_cat = cat(3,thisBA6(:,1:20),thisBA6(:,21:40),thisBA6(:,41:60),thisBA6(:,61:80));
elseif gridX == 4
    thisBA6_cat = cat(3,thisBA6(:,1:4),thisBA6(:,5:8),thisBA6(:,9:12),thisBA6(:,13:16));
end
thisBA6_cat_mean = mean(thisBA6_cat,3);
% now smooth
%[thisBA2_cat_mean_smoothed] = smoothMyRF(thisBA2_cat_mean);

figure('Position',[100 100 1600 400])
tiledlayout(2,4)
nexttile
imagesc(thisBA3a_cat_mean)
[x0, y0, v] = find2dMax(thisBA3a_cat_mean);
hold on
s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
    'markeredgecolor', 'k');
t_ = text(10,-3,'BA3a');
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
t_ = text(10,-3,'BA3b');
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
t_ = text(10,-3,'BA1');
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
t_ = text(10,-3,'BA2');
set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
    'helvetica', 'interpreter', 'none','verticalalignment','top');
ax = gca;
ax.FontSize = 12;
colorbar
axis square
colormap jet

nexttile
imagesc(thisBA4a_cat_mean)
[x0, y0, v] = find2dMax(thisBA4a_cat_mean);
hold on
s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
    'markeredgecolor', 'k');
t_ = text(10,-3,'BA4a');
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
t_ = text(10,-3,'BA4p');
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
t_ = text(10,-3,'BA6');
set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
    'helvetica', 'interpreter', 'none','verticalalignment','top');
ax = gca;
ax.FontSize = 12;
colorbar
axis square
colormap jet


print(gcf, sprintf('grandmeanofPCA_summed_weighted_%s_PC%s_BAA_interped_BIG',currentModel, num2str(noComponents)), '-r300', '-dpng')

% % INTERP
% figure('Position',[100 100 1600 400])
% tiledlayout(1,4)
% nexttile
% imagesc(thisBA3a_cat_mean_smoothed)
% [x0, y0, v] = find2dMax(thisBA3a_cat_mean);
% hold on
% s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%     'markeredgecolor', 'k');
% t_ = text(10,-1,'BA3a');
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
% t_ = text(10,-1,'BA3b');
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
% t_ = text(10,-1,'BA1');
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
% t_ = text(10,-1,'BA2');
% set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
%     'helvetica', 'interpreter', 'none','verticalalignment','top');
% ax = gca;
% ax.FontSize = 12;
% colorbar
% axis square
%
%
% print(gcf, sprintf('grandmeanofPCA_summed_weighted_%s_PC%s_BAA_interp',currentModel, num2str(noComponents)), '-r300', '-dpng')
%
% %
%% BASETIPS VERSION
%
%
for iSub = 1:8
    for iBA = 1:4
        for iDig = 1:4
            
            %size error check
            if size(theGRANDpca_basetips.myGRANDpca.group{iSub}{iBA,iDig},2) < length(noComponents)
                thisComponent = 1;
            else
                thisComponent = noComponents;
            end
            
            digits.location{iBA,iDig,iSub} = theGRANDpca_basetips.myGRANDpca.group{iSub}{iBA,iDig}(:,thisComponent);
            digits.explained{iBA,iDig,iSub} = theGRANDpca_basetips.myGRANDpca.explained{iSub}{iBA,iDig}(thisComponent,:);
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
    for iBA = 1:4
        for iDig = 1:4
            
            %size error check
            if size(digits.location{iBA,iDig,iSub},2) < length(noComponentidx)
                thisComponent = 1;
            else
                thisComponent = noComponentidx;
            end
            
            digits.weighted{iBA,iDig,iSub} = digits.location{iBA,iDig,iSub}(:,thisComponent) * (digits.explained{iBA,iDig,iSub}(thisComponent)./100);
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
        digits.meanPCAunweightedRFreshaped{iBA,iDig} = reshape(digits.meanPCAunweightedRF{iBA,iDig},gridX,gridY);
    end
end

theDigs = 2:5;

% now display
figure('Position',[100 100 1200 1000])
tiledlayout(4,4)
theBaa = {'PD2','PD3','PD4','PD5'};
for iBa = 1:4
    myBa = theBaa{iBa};
    
    for iDig = 1:4
        myDig = theDigs(iDig);
        nexttile
        imagesc(digits.meanPCAweightedRFreshaped{iBa,iDig});
        [x0, y0, v] = find2dMax(digits.meanPCAweightedRFreshaped{iBa,iDig});
        hold on
        s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
            'markeredgecolor', 'k');
        t_ = text(0,0,sprintf('dig:%d, %s, model:%s', myDig, ...
            myBa, currentModel));
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

print(gcf, sprintf('grandmeanofPCA_weighted_%s_PC%s_BASETIPS',currentModel, num2str(noComponents)), '-r300', '-dpng')

%% try summing
thisgrid = digits.meanPCAweightedRFreshaped;
thisBA3a = cell2mat(thisgrid(1,:)); %D2
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
% now smooth
%[thisBA3b_cat_mean_smoothed] = smoothMyRF(thisBA3b_cat_mean);

thisBA1 = cell2mat(thisgrid(3,:));%D4
if gridX == 20
    thisBA1_cat = cat(3,thisBA1(:,1:20),thisBA1(:,21:40),thisBA1(:,41:60),thisBA1(:,61:80));
elseif gridX == 4
    thisBA1_cat = cat(3,thisBA1(:,1:4),thisBA1(:,5:8),thisBA1(:,9:12),thisBA1(:,13:16));
end
thisBA1_cat_mean = mean(thisBA1_cat,3);
% now smooth
%[thisBA1_cat_mean_smoothed] = smoothMyRF(thisBA1_cat_mean);

thisBA2 = cell2mat(thisgrid(4,:));%D5
if gridX == 20
    thisBA2_cat = cat(3,thisBA2(:,1:20),thisBA2(:,21:40),thisBA2(:,41:60),thisBA2(:,61:80));
elseif gridX == 4
    thisBA2_cat = cat(3,thisBA2(:,1:4),thisBA2(:,5:8),thisBA2(:,9:12),thisBA2(:,13:16));
end
thisBA2_cat_mean = mean(thisBA2_cat,3);
% now smooth
%[thisBA2_cat_mean_smoothed] = smoothMyRF(thisBA2_cat_mean);

% raw version
figure('Position',[100 100 1600 400])
tiledlayout(1,4)
nexttile
imagesc(thisBA3a_cat_mean)
[x0, y0, v] = find2dMax(thisBA3a_cat_mean);
hold on
s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
    'markeredgecolor', 'k');
t_ = text(10,-1,'PD1');
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
t_ = text(10,-1,'PD2');
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
t_ = text(10,-1,'PD3');
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
t_ = text(10,-1,'PD4');
set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
    'helvetica', 'interpreter', 'none','verticalalignment','top');
ax = gca;
ax.FontSize = 12;
colorbar
axis square
colormap jet

print(gcf, sprintf('grandmeanofPCA_summed_weighted_%s_PC%s_BASETIPS_interped',currentModel, num2str(noComponents)), '-r300', '-dpng')

% % INTERP VERSION
%
% figure('Position',[100 100 1600 400])
% tiledlayout(1,4)
% nexttile
% imagesc(thisBA3a_cat_mean_smoothed)
% [x0, y0, v] = find2dMax(thisBA3a_cat_mean);
% hold on
% s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%     'markeredgecolor', 'k');
% t_ = text(10,-1,'PD1');
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
% t_ = text(10,-1,'PD2');
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
% t_ = text(10,-1,'PD3');
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
% t_ = text(10,-1,'PD4');
% set(t_, 'color', [0 0 0 ],'fontsize', 16, 'fontname', ...
%     'helvetica', 'interpreter', 'none','verticalalignment','top');
% ax = gca;
% ax.FontSize = 12;
% colorbar
% axis square
%
% print(gcf, sprintf('grandmeanofPCA_summed_weighted_%s_PC%s_BASETIPS_interp',currentModel, num2str(noComponents)), '-r300', '-dpng')
%
%


%%

% UNWEIGHTED VERSION

% figure('Position',[100 100 1000 1000])
% tiledlayout(4,4)
% for iBa = 1:4
%     for iDig = 1:4
%         nexttile
%         imagesc(digits.meanPCAunweightedRFreshaped{iBa,iDig});
%         [x0, y0, v] = find2dMax(digits.meanPCAunweightedRFreshaped{iBa,iDig});
%         hold on
%         s_= scatter(x0,y0,'filled','markerfacecolor', 'w', ...
%             'markeredgecolor', 'k');
%         t_ = text(0,0,sprintf('dig:%d, pd%d, model:%s', iDig, ...
%             iBa, currentModel));
%         set(t_, 'color', [0 0 0 ], ...
%             'fontsize', 12, 'fontname', 'helvetica', 'interpreter', 'none', ...
%             'verticalalignment', 'top');
%         axis square
%     end
% end
%
% print(gcf, sprintf('grandmeanofPCA_unweighted_%s_PC%s',currentModel, num2str(noComponents)), '-r300', '-dpng')
%
%



%% Option C
% run PCA directly on stacked all subjects RFs
%
% mypath = '/Volumes/nemosine/prfsomato_fits';
% cd(mypath)
% close all
% clc
%
% a = 1; % principal component
% %noComponents = 16;
%
% theGRANDrfs = load('prf_2ddouble_dd_nov_myRFs.mat');
%
% for iSub = 1:10
%
%     d2pd2.group{iSub} = theGRANDrfs.myRF.group{iSub}{1,1};
%     d2pd3.group{iSub} = theGRANDrfs.myRF.group{iSub}{1,2};
%     d2pd4.group{iSub} = theGRANDrfs.myRF.group{iSub}{1,3};
%     d2pd5.group{iSub} = theGRANDrfs.myRF.group{iSub}{1,4};
%
%     d3pd2.group{iSub} = theGRANDrfs.myRF.group{iSub}{2,1};
%     d3pd3.group{iSub} = theGRANDrfs.myRF.group{iSub}{2,2};
%     d3pd4.group{iSub} = theGRANDrfs.myRF.group{iSub}{2,3};
%     d3pd5.group{iSub} = theGRANDrfs.myRF.group{iSub}{2,4};
%
%     d4pd2.group{iSub} = theGRANDrfs.myRF.group{iSub}{3,1};
%     d4pd3.group{iSub} = theGRANDrfs.myRF.group{iSub}{3,2};
%     d4pd4.group{iSub} = theGRANDrfs.myRF.group{iSub}{3,3};
%     d4pd5.group{iSub} = theGRANDrfs.myRF.group{iSub}{3,4};
%
%     d5pd2.group{iSub} = theGRANDrfs.myRF.group{iSub}{4,1};
%     d5pd3.group{iSub} = theGRANDrfs.myRF.group{iSub}{4,2};
%     d5pd4.group{iSub} = theGRANDrfs.myRF.group{iSub}{4,3};
%     d5pd5.group{iSub} = theGRANDrfs.myRF.group{iSub}{4,4};
%
% end
%
% d2pd2_bigRF = cat(3,d2pd2.group{1},d2pd2.group{2},d2pd2.group{3},d2pd2.group{4},d2pd2.group{5},d2pd2.group{6},d2pd2.group{7},d2pd2.group{8},d2pd2.group{9},d2pd2.group{10});
% d2pd3_bigRF = cat(3,d2pd3.group{1},d2pd3.group{2},d2pd3.group{3},d2pd3.group{4},d2pd3.group{5},d2pd3.group{6},d2pd3.group{7},d2pd3.group{8},d2pd3.group{9},d2pd3.group{10});
% d2pd4_bigRF = cat(3,d2pd4.group{1},d2pd4.group{2},d2pd4.group{3},d2pd4.group{4},d2pd4.group{5},d2pd4.group{6},d2pd4.group{7},d2pd4.group{8},d2pd4.group{9},d2pd4.group{10});
% d2pd5_bigRF = cat(3,d2pd5.group{1},d2pd5.group{2},d2pd5.group{3},d2pd5.group{4},d2pd5.group{5},d2pd5.group{6},d2pd5.group{7},d2pd5.group{8},d2pd5.group{9},d2pd5.group{10});
% d3pd2_bigRF = cat(3,d3pd2.group{1},d3pd2.group{2},d3pd2.group{3},d3pd2.group{4},d2pd2.group{5},d3pd2.group{6},d3pd2.group{7},d3pd2.group{8},d3pd2.group{9},d3pd2.group{10});
% d3pd3_bigRF = cat(3,d3pd3.group{1},d3pd3.group{2},d3pd3.group{3},d3pd3.group{4},d2pd3.group{5},d3pd3.group{6},d3pd3.group{7},d3pd3.group{8},d3pd3.group{9},d3pd3.group{10});
% d3pd4_bigRF = cat(3,d3pd4.group{1},d3pd4.group{2},d3pd4.group{3},d3pd4.group{4},d2pd4.group{5},d3pd4.group{6},d3pd4.group{7},d3pd4.group{8},d3pd4.group{9},d3pd4.group{10});
% d3pd5_bigRF = cat(3,d3pd5.group{1},d3pd5.group{2},d3pd5.group{3},d3pd5.group{4},d2pd5.group{5},d3pd5.group{6},d3pd5.group{7},d3pd5.group{8},d3pd5.group{9},d3pd5.group{10});
% d4pd2_bigRF = cat(3,d4pd2.group{1},d4pd2.group{2},d4pd2.group{3},d4pd2.group{4},d2pd2.group{5},d4pd2.group{6},d4pd2.group{7},d4pd2.group{8},d4pd2.group{9},d4pd2.group{10});
% d4pd3_bigRF = cat(3,d4pd3.group{1},d4pd3.group{2},d4pd3.group{3},d4pd3.group{4},d2pd3.group{5},d4pd3.group{6},d4pd3.group{7},d4pd3.group{8},d4pd3.group{9},d4pd3.group{10});
% d4pd4_bigRF = cat(3,d4pd4.group{1},d4pd4.group{2},d4pd4.group{3},d4pd4.group{4},d2pd4.group{5},d4pd4.group{6},d4pd4.group{7},d4pd4.group{8},d4pd4.group{9},d4pd4.group{10});
% d4pd5_bigRF = cat(3,d4pd5.group{1},d4pd5.group{2},d4pd5.group{3},d4pd5.group{4},d2pd5.group{5},d4pd5.group{6},d4pd5.group{7},d4pd5.group{8},d4pd5.group{9},d4pd5.group{10});
% d5pd2_bigRF = cat(3,d5pd2.group{1},d5pd2.group{2},d5pd2.group{3},d5pd2.group{4},d2pd2.group{5},d5pd2.group{6},d5pd2.group{7},d5pd2.group{8},d5pd2.group{9},d5pd2.group{10});
% d5pd3_bigRF = cat(3,d5pd3.group{1},d5pd3.group{2},d5pd3.group{3},d5pd3.group{4},d2pd3.group{5},d5pd3.group{6},d5pd3.group{7},d5pd3.group{8},d5pd3.group{9},d5pd3.group{10});
% d5pd4_bigRF = cat(3,d5pd4.group{1},d5pd4.group{2},d5pd4.group{3},d5pd4.group{4},d2pd4.group{5},d5pd4.group{6},d5pd4.group{7},d5pd4.group{8},d5pd4.group{9},d5pd4.group{10});
% d5pd5_bigRF = cat(3,d5pd5.group{1},d5pd5.group{2},d5pd5.group{3},d5pd5.group{4},d2pd5.group{5},d5pd5.group{6},d5pd5.group{7},d5pd5.group{8},d5pd5.group{9},d5pd5.group{10});
%
% % now PCA these
% d2pd2_X = reshape(d2pd2_bigRF, [size(d2pd2_bigRF,1).*size(d2pd2_bigRF,1) size(d2pd2_bigRF,3)]);
% % need to flip it as ROWS = Observations, and COLUMNS = Variables;
% % so e.g. we have 48 observations (voxels), and 16 variables (each pRF
% % shape), for one ROI, e.g. D2, PD2
% d2pd3_X = reshape(d2pd3_bigRF, [size(d2pd3_bigRF,1).*size(d2pd3_bigRF,1) size(d2pd3_bigRF,3)]);
% d2pd4_X = reshape(d2pd4_bigRF, [size(d2pd4_bigRF,1).*size(d2pd4_bigRF,1) size(d2pd4_bigRF,3)]);
% d2pd5_X = reshape(d2pd5_bigRF, [size(d2pd5_bigRF,1).*size(d2pd5_bigRF,1) size(d2pd5_bigRF,3)]);
%
% d3pd2_X = reshape(d3pd2_bigRF, [size(d3pd2_bigRF,1).*size(d3pd2_bigRF,1) size(d3pd2_bigRF,3)]);
% d3pd3_X = reshape(d3pd3_bigRF, [size(d3pd3_bigRF,1).*size(d3pd3_bigRF,1) size(d3pd3_bigRF,3)]);
% d3pd4_X = reshape(d3pd4_bigRF, [size(d3pd4_bigRF,1).*size(d3pd4_bigRF,1) size(d3pd4_bigRF,3)]);
% d3pd5_X = reshape(d3pd5_bigRF, [size(d3pd5_bigRF,1).*size(d3pd5_bigRF,1) size(d3pd5_bigRF,3)]);
%
% d4pd2_X = reshape(d4pd2_bigRF, [size(d4pd2_bigRF,1).*size(d4pd2_bigRF,1) size(d4pd2_bigRF,3)]);
% d4pd3_X = reshape(d4pd3_bigRF, [size(d4pd3_bigRF,1).*size(d4pd3_bigRF,1) size(d4pd3_bigRF,3)]);
% d4pd4_X = reshape(d4pd4_bigRF, [size(d4pd4_bigRF,1).*size(d4pd4_bigRF,1) size(d4pd4_bigRF,3)]);
% d4pd5_X = reshape(d4pd5_bigRF, [size(d4pd5_bigRF,1).*size(d4pd5_bigRF,1) size(d4pd5_bigRF,3)]);
%
% d5pd2_X = reshape(d5pd2_bigRF, [size(d5pd2_bigRF,1).*size(d5pd2_bigRF,1) size(d5pd2_bigRF,3)]);
% d5pd3_X = reshape(d5pd3_bigRF, [size(d5pd3_bigRF,1).*size(d5pd3_bigRF,1) size(d5pd3_bigRF,3)]);
% d5pd4_X = reshape(d5pd4_bigRF, [size(d5pd4_bigRF,1).*size(d5pd4_bigRF,1) size(d5pd4_bigRF,3)]);
% d5pd5_X = reshape(d5pd5_bigRF, [size(d5pd5_bigRF,1).*size(d5pd5_bigRF,1) size(d5pd5_bigRF,3)]);
%
% % transpose
% Xt_d2pd2 = transpose(d2pd2_X);
% Xt_d2pd3 = transpose(d2pd3_X);
% Xt_d2pd4 = transpose(d2pd4_X);
% Xt_d2pd5 = transpose(d2pd5_X);
%
% Xt_d3pd2 = transpose(d3pd2_X);
% Xt_d3pd3 = transpose(d3pd3_X);
% Xt_d3pd4 = transpose(d3pd4_X);
% Xt_d3pd5 = transpose(d3pd5_X);
%
% Xt_d4pd2 = transpose(d4pd2_X);
% Xt_d4pd3 = transpose(d4pd3_X);
% Xt_d4pd4 = transpose(d4pd4_X);
% Xt_d4pd5 = transpose(d4pd5_X);
%
% Xt_d5pd2 = transpose(d5pd2_X);
% Xt_d5pd3 = transpose(d5pd3_X);
% Xt_d5pd4 = transpose(d5pd4_X);
% Xt_d5pd5 = transpose(d5pd5_X);
%
%
%
%
%
% % PCA
% [coeff_d2pd2] = pca(Xt_d2pd2);
% [coeff_d2pd3] = pca(Xt_d2pd3);
% [coeff_d2pd4] = pca(Xt_d2pd4);
% [coeff_d2pd5] = pca(Xt_d2pd5);
%
% [coeff_d3pd2] = pca(Xt_d3pd2);
% [coeff_d3pd3] = pca(Xt_d3pd3);
% [coeff_d3pd4] = pca(Xt_d3pd4);
% [coeff_d3pd5] = pca(Xt_d3pd5);
%
% [coeff_d4pd2] = pca(Xt_d4pd2);
% [coeff_d4pd3] = pca(Xt_d4pd3);
% [coeff_d4pd4] = pca(Xt_d4pd4);
% [coeff_d4pd5] = pca(Xt_d4pd5);
%
% [coeff_d5pd2] = pca(Xt_d5pd2);
% [coeff_d5pd3] = pca(Xt_d5pd3);
% [coeff_d5pd4] = pca(Xt_d5pd4);
% [coeff_d5pd5] = pca(Xt_d5pd5);
%
% % reshape and plot
% figure('Position',[100 100 1000 1000])
%
% tiledlayout(4,4)
% nexttile
% imagesc(reshape(coeff_d2pd2(:,a),4,4))
% title('D2PD2')
% axis square
% nexttile
% imagesc(reshape(coeff_d3pd2(:,a),4,4))
% title('D3PD2')
% axis square
% nexttile
% imagesc(reshape(coeff_d4pd2(:,a),4,4))
% title('D4PD2')
% axis square
% nexttile
% imagesc(reshape(coeff_d5pd2(:,a),4,4))
% title('D5PD2')
% axis square
%
% nexttile
% imagesc(reshape(coeff_d2pd3(:,a),4,4))
% title('D2PD3')
% axis square
% nexttile
% imagesc(reshape(coeff_d3pd3(:,a),4,4))
% title('D3PD3')
% axis square
% nexttile
% imagesc(reshape(coeff_d4pd3(:,a),4,4))
% title('D4PD3')
% axis square
% nexttile
% imagesc(reshape(coeff_d5pd3(:,a),4,4))
% title('D5PD3')
% axis square
%
% nexttile
% imagesc(reshape(coeff_d2pd4(:,a),4,4))
% title('D2PD4')
% axis square
% nexttile
% imagesc(reshape(coeff_d3pd4(:,a),4,4))
% title('D3PD4')
% axis square
% nexttile
% imagesc(reshape(coeff_d4pd4(:,a),4,4))
% title('D4PD4')
% axis square
% nexttile
% imagesc(reshape(coeff_d5pd4(:,a),4,4))
% title('D5PD4')
% axis square
%
% nexttile
% imagesc(reshape(coeff_d2pd5(:,a),4,4))
% title('D2PD5')
% axis square
% nexttile
% imagesc(reshape(coeff_d3pd5(:,a),4,4))
% title('D3PD5')
% axis square
% nexttile
% imagesc(reshape(coeff_d4pd5(:,a),4,4))
% title('D4PD5')
% axis square
% nexttile
% imagesc(reshape(coeff_d5pd5(:,a),4,4))
% title('D5PD5')
% axis square
%
% print(gcf, sprintf('grandPCA_fromRFs'), '-r300', '-dpng')



%% do means
mypath='pRF_somato_data/data/';
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




