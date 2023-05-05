clc
close all


mypath = '/Volumes/nemosine/prfsomato_fits';
cd(mypath)


mysubs = {'prf1','prf2','prf3','prf4','prf6','prf7','prf8','prf10'};
digits = 2:5;
locs = 2:8;

datafilenames = {'prf_2ddouble_dd_nov_big', ...
    'prf_1d_nodiag2_nov_big', ...
    'prf_1dt_nodiag2_nov_big', ...
    'prf_sixteen_nodiag2_com_nov_big'};

thisModel = 1;
%currentModel = datafilenames{thisModel};

models = {'2DGauss','1DWD-Gauss','1DBD-Gauss','Uncon'};
%mymodel = models{thisModel};

if thisModel == 1
    load('prf_2ddouble_dd_nov_bigBAA_myBestTC.mat');
elseif thisModel == 2
    load('prf_1d_nodiag2_nov_bigBAA_myBestTC.mat');
elseif thisModel == 3
    load('prf_1dt_nodiag2_nov_bigBAA_myBestTC.mat');
elseif thisModel == 4
    load('prf_sixteen_nodiag2_com_nov_bigBAA_myBestTC.mat');
end

% models = {'2DGauss','1DWD-Gauss','1DBD-Gauss','Uncon'};
% mymodel = models{thisModel};
thedigits = {'D2','D3','D4','D5'};
%mydigit = thedigits{thisDig};
areas = {'BA3a','BA3b','BA1','BA2','BA4a','BA4p','BA6'};
%myarea = areas{thisLoc};


%
tic
for iSub = 1:length(mysubs)
    %figure('Position',[100 100 1200 1200])
    %tiledlayout(7,4)
    currentSub = mysubs{iSub};
    currentModel = models{thisModel};

    for iLoc = 1:7
        %currentLoc = locs(iLoc);
        currentLoc = areas{iLoc};

        for iDigit = 1:4

            %currentDigit = digits(iDigit);
            currentDigit = thedigits{iDigit};

            %nexttile
            
            segmentwh = myBestTC.tc{iSub}{iLoc,iDigit};
            segment = segmentwh(169:169+95,:);
            
            segmentmo = myBestTC.modresp{iSub}{iLoc,iDigit};
            segmo = segmentmo(169:169+95,:);
            
            themeanR2(iLoc,iDigit,iSub) = mean(myBestTC.bestr2{iSub}{iLoc,iDigit});
            
            choppedsegment = cat(3,segment(1:8,:),segment(9:16,:),segment(17:24,:),segment(25:32,:),segment(33:40,:),...
                segment(41:48,:),segment(49:56,:),segment(57:64,:),segment(65:72,:),segment(73:80,:),...
                segment(81:88,:),segment(89:96,:));
            
            choppedsegmo = cat(3,segmo(1:8,:),segmo(9:16,:),segmo(17:24,:),segmo(25:32,:),segmo(33:40,:),...
                segmo(41:48,:),segmo(49:56,:),segmo(57:64,:),segmo(65:72,:),segmo(73:80,:),...
                segmo(81:88,:),segmo(89:96,:));

            
            choppedsegmentM = mean(choppedsegment,3);
            choppedsegmentMM(:,iLoc,iDigit,iSub) = mean(choppedsegmentM,2);
            
            choppedsegmentSE = std(choppedsegmentM,[],2)./sqrt(length(choppedsegmentM));
            choppedsegmentSE(:,iLoc,iDigit,iSub) = choppedsegmentSE(:);
            choppedsegmoM = mean(choppedsegmo,3);
            choppedsegmoMM(:,iLoc,iDigit,iSub) = mean(choppedsegmoM,2);


            %figure('Position',[100 100 1200 600],'PaperOrientation','landscape')
            %tiledlayout(2,3)
            %nexttile([1 2])
            %figure('Position',[100 100 600 300])
%             plot(segment,'linewidth',1.5,'color','k')
%             hold on
%             plot(segmo,'linewidth',0.8,'color','r','Marker','o')
%             %ylim([0.92 1.08])
%             %text(10,0.95,['Adj r2=' num2str(myBestTC.bestr2{thisSub}{thisLoc,thisDig}(:,whichvox))]);
%             xlabel('Time (s)')
%             ylabel('Normalised Bold %')
%             %title(sprintf('Good fitting voxel %s %s %s idx%d',mydigit, myarea, mymodel,thisGuy))
%             

% ------------------------------------------------------------------------------------

           % this is the old pltoting for one subject 
%             errorbar(choppedsegmentMM,choppedsegmentSE,'linewidth',1.5,'color','k')
%             hold on
%             plot(choppedsegmoMM,'linewidth',1.5,'color','r','Marker','o')
%             %title('Mean & SE','FontSize',10)
%             title(sprintf('sub:%s, dig:%s, pd%s, m:%s ', currentSub, currentDigit, ...
%                  currentLoc, currentModel), 'FontSize',10)
%             xlabel('Time (s)','FontSize',10)
%             ylabel('Normalised Bold %','FontSize',10)
%     
%             text(0,1,['Adj r2=' num2str(themeanR2)],'fontsize', 10);
% 
% 
%             a = get(gca,'XTickLabel');
%             set(gca,'XTickLabel',a,'fontsize',10)
%             set(gca,'XTickLabelMode','auto')
%             
 

% ------------------------------------------------------------------------------------
            
            
            

            %xlim([1 8])
%             t_ = text(1,1,sprintf('sub:%s, dig:%d, pd%d, m:%s', currentSub, currentDigit, ...
%                  currentLoc, currentModel));
%             set(t_, 'color', [0 0 0 ], ...
%                 'fontsize', 10, 'fontname', 'helvetica', 'interpreter', 'none', ...
%                 'verticalalignment', 'top');



            %plot(myBestTC.tc{iSub}{iLoc,iDigit},'linewidth',1.5)
            %hold on
            %plot(myBestTC.modresp{iSub}{iLoc,iDigit})

%            ylim([0.92 1.08])

%             t_ = text(0,1.1,sprintf('sub:%s, dig:%d, pd%d, m:%s', currentSub, currentDigit, ...
%                 currentLoc, currentModel));
%             set(t_, 'color', [0 0 0 ], ...
%                 'fontsize', 10, 'fontname', 'helvetica', 'interpreter', 'none', ...
%                 'verticalalignment', 'top');
%             hold off

        end
    end
end
toc
% now mean across subjects

% choppedsegmentMM is the TC, size = 8 7 4 8 (tc, locs, digs, subs)

subCHOP_tc = mean(choppedsegmentMM,4);
subCHOP_se = mean(choppedsegmentSE,4);
subCHOP_mo = mean(choppedsegmoMM,4);
subCHOP_r2 = mean(themeanR2,3);
% these r2 actually represent the full timecourse, not a single cycle

% r2 for single cycle, 1 model
for iLoc = 1:7
    for iDig = 1:4

        SSR(iLoc,iDig) = sum((subCHOP_mo(:,iLoc,iDig) - subCHOP_tc(:,iLoc,iDig)).^2);
        TSS(iLoc,iDig) = sum((subCHOP_tc(:,iLoc,iDig) - mean(subCHOP_tc(:,iLoc,iDig))).^2);
        Rsquared(iLoc,iDig) = 1-SSR(iLoc,iDig)/TSS(iLoc,iDig);

    end
end

%
figure('Position',[100 100 600 1200])
tiledlayout(7,4)

for iLoc = 1:7
    currentLoc = areas{iLoc};
    for iDig = 1:4
        currentDigit = thedigits{iDig};
        nexttile
        errorbar(subCHOP_tc(:,iLoc,iDig),subCHOP_se(:,iLoc,iDig),'linewidth',1.5,'color','k')
        xlim([1 8])
        hold on
        plot(subCHOP_mo(:,iLoc,iDig),'linewidth',1.5,'color','r','Marker','o')
        %title(sprintf('dig:%s, pd%s, m:%s', currentDigit, ...
        %    currentLoc, currentModel), 'FontSize',10)
        %xlabel('Time (TRs)','FontSize',10)
        %ylabel('Normalised Bold %','FontSize',10)
        text(1.5,1.008,sprintf('%.2g',Rsquared(iLoc,iDig)),'fontsize', 10);
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',10)
        set(gca,'XTickLabelMode','auto')
        axis tight
        ylim([0.99 1.01])
        yticks([0.99:0.005:1.01])
        grid MINOR
    end
end



%
mymodel = models{thisModel};
filename =  sprintf('TC_%s.pdf',mymodel);
filename2 = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',filename);
%export_fig(filename2,'-pdf');
print(gcf,filename2, '-r300', '-dpdf','-bestfit')

%% now plot the means collapsed across digits

msubCHOP_tc = mean(subCHOP_tc,3);
msubCHOP_se = std(subCHOP_tc,[],3);
msubCHOP_mo = mean(subCHOP_mo,3);


figure('Position',[100 100 150 1200])
tiledlayout(7,1)
for ii = 1:7
    nexttile
    %errorbar(msubCHOP_tc(:,ii),msubCHOP_se(:,ii),'linewidth',1.5,'color','k')
    shadedErrorBar([],msubCHOP_tc(:,ii),msubCHOP_se(:,ii),'lineProps','-k','transparent',1);
    xlim([1 8])
    hold on
    plot(msubCHOP_mo(:,ii),'linewidth',1,'color','r','Marker','o')

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',10)
    set(gca,'XTickLabelMode','auto')
    axis tight
    ylim([0.99 1.01])
    %yticks([0.99:0.005:1.01])
    xticklabels({})
    yticklabels({})
    grid MINOR
end
filenameX =  sprintf('TC_%s_average.pdf',mymodel);
filenameXX = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',filenameX);
%export_fig(filename2,'-pdf');
print(gcf,filenameXX, '-r300', '-dpdf','-bestfit')


%% mean across digits
% 
% figure
% tiledlayout(7,1)
% for ii = 1:7
%     cattest = cat(2,myBestTC.tc{1}{ii,1},myBestTC.tc{1}{ii,2},myBestTC.tc{1}{ii,3},myBestTC.tc{1}{ii,4});
%     catmean(:,ii) = mean(cattest,2);
%     nexttile
%     plot(catmean(:,ii))
% end




%%
% for iSub = 1
%     figure('Position',[100 100 1800 1200])
%     tiledlayout(7,4)
%     currentSub = mysubs{iSub};
%     for iLoc = 1:7
%         currentLoc = locs(iLoc);
%         for iDigit = 1:4
%             
%             currentDigit = digits(iDigit);
%             nexttile
%             
%             plot(myBestTC.tc{iSub}{iLoc,iDigit},'linewidth',1.5)
%             %hold on
%             %plot(myBestTC.modresp{iSub}{iLoc,iDigit})
%             
%             ylim([0.92 1.08])
%             
%             t_ = text(0,1.1,sprintf('sub:%s, dig:%d, pd%d, m:%s', currentSub, currentDigit, ...
%                 currentLoc, currentModel));
%             set(t_, 'color', [0 0 0 ], ...
%                 'fontsize', 10, 'fontname', 'helvetica', 'interpreter', 'none', ...
%                 'verticalalignment', 'top');
%             hold off
%             
%         end
%     end
% end
% 
% %% mean across digits
% 
% figure
% tiledlayout(7,1)
% for ii = 1:7
%     cattest = cat(2,myBestTC.tc{1}{ii,1},myBestTC.tc{1}{ii,2},myBestTC.tc{1}{ii,3},myBestTC.tc{1}{ii,4});
%     catmean(:,ii) = mean(cattest,2);
%     nexttile
%     plot(catmean(:,ii))
% end


%
%close all 
thisLoc = 3;
thisDig = 3;

thedigits = {'D2','D3','D4','D5'};
mydigit = thedigits{thisDig};
areas = {'BA3a','BA3b','BA1','BA2'};
myarea = areas{thisLoc};

thisLocB = 4;
thisDigB = 4;

mydigitB = thedigits{thisDigB};

myareaB = areas{thisLocB};


thisSub = 1;

thisGuy = 1233;
idxvox = ismember(myBestTC.MLI{thisSub}{thisLoc,thisDig},thisGuy);
whichvox = find(idxvox);

thisGuyB = 1884;
idxvoxB = ismember(myBestTC.MLI{thisSub}{thisLocB,thisDigB},thisGuyB);
whichvoxB = find(idxvoxB);

%whichvox = 4;
% whichvoxB = 10;

%myidx = num2str(myBestTC.MLI{thisSub}{thisLoc,thisDig}(whichvox));

models = {'2DGauss','1DWD-Gauss','1DBD-Gauss','Uncon'};
mymodel = models{thisModel};


segmentwh = myBestTC.tc{thisSub}{thisLoc,thisDig}(:,whichvox);
segment = segmentwh(169:169+95);

segmentmo = myBestTC.modresp{thisSub}{thisLoc,thisDig}(:,whichvox);
segmo = segmentmo(169:169+95);

segmentwhbad = myBestTC.tc{thisSub}{thisLocB,thisDigB}(:,whichvoxB);
segmentbad = segmentwhbad(169:169+95);

segmentmobad = myBestTC.modresp{thisSub}{thisLocB,thisDigB}(:,whichvoxB);
segmobad = segmentmobad(169:169+95);

%try chop
% 96 / 12 = 8

% good
choppedsegment = cat(2,segment(1:8),segment(9:16),segment(17:24),segment(25:32),segment(33:40),...
    segment(41:48),segment(49:56),segment(57:64),segment(65:72),segment(73:80),...
    segment(81:88),segment(89:96));

choppedsegmo = cat(2,segmo(1:8),segmo(9:16),segmo(17:24),segmo(25:32),segmo(33:40),...
    segmo(41:48),segmo(49:56),segmo(57:64),segmo(65:72),segmo(73:80),...
    segmo(81:88),segmo(89:96));

choppedsegmentM = mean(choppedsegment,2);
choppedsegmentSE = std(choppedsegment,[],2)./sqrt(length(choppedsegment));
choppedsegmentSE = choppedsegmentSE(:);
choppedsegmoM = mean(choppedsegmo,2);

% bad
choppedsegmentB = cat(2,segmentbad(1:8),segmentbad(9:16),segmentbad(17:24),segmentbad(25:32),segmentbad(33:40),...
    segmentbad(41:48),segmentbad(49:56),segmentbad(57:64),segmentbad(65:72),segmentbad(73:80),...
    segmentbad(81:88),segmentbad(89:96));

choppedsegmoB = cat(2,segmobad(1:8),segmobad(9:16),segmobad(17:24),segmobad(25:32),segmobad(33:40),...
    segmobad(41:48),segmobad(49:56),segmobad(57:64),segmobad(65:72),segmobad(73:80),...
    segmobad(81:88),segmobad(89:96));

choppedsegmentMB = mean(choppedsegmentB,2);
choppedsegmentSEB = std(choppedsegmentB,[],2)./sqrt(length(choppedsegmentB));
choppedsegmentSEB = choppedsegmentSEB(:);
choppedsegmoMB = mean(choppedsegmoB,2);


%
figure('Position',[100 100 1200 600],'PaperOrientation','landscape')
tiledlayout(2,3)
nexttile([1 2])
%figure('Position',[100 100 600 300])
plot(segment,'linewidth',1.5,'color','k')
hold on
plot(segmo,'linewidth',0.8,'color','r','Marker','o')
ylim([0.92 1.08])
text(10,0.95,['Adj r2=' num2str(myBestTC.bestr2{thisSub}{thisLoc,thisDig}(:,whichvox))]);
xlabel('Time (s)')
ylabel('Normalised Bold %')
title(sprintf('Good fitting voxel %s %s %s idx%d',mydigit, myarea, mymodel,thisGuy))

nexttile
errorbar(choppedsegmentM,choppedsegmentSE,'linewidth',1.5,'color','k')
hold on
plot(choppedsegmoM,'linewidth',1.5,'color','r','Marker','o')
title('Mean & SE')
xlabel('Time (s)')
ylabel('Normalised Bold %')
xlim([1 8])

nexttile([1 2])
plot(segmentbad,'linewidth',1.5,'color','k')
hold on
plot(segmobad,'linewidth',0.8,'color','r','Marker','o')
ylim([0.92 1.08])
text(10,0.95,['Adj r2=' num2str(myBestTC.bestr2{thisSub}{thisLocB,thisDigB}(:,whichvoxB))]);
xlabel('Time (s)')
ylabel('Normalised Bold %')
title(sprintf('Poor fitting voxel %s %s %s idx%d',mydigitB, myareaB, mymodel,thisGuyB))

nexttile
errorbar(choppedsegmentMB,choppedsegmentSEB,'linewidth',1.5,'color','k')
hold on
plot(choppedsegmoMB,'linewidth',1.5,'color','r','Marker','o')
title('Mean & SE')
xlabel('Time (s)')
ylabel('Normalised Bold %')
xlim([1 8])

print(gcf,sprintf('TC_%s.pdf',mymodel), '-r300', '-dpdf','-bestfit')



% 
% figure
% plot(myBestTC.tc{1}{3,3}(:,1),'linewidth',1.5)
% hold on
% plot(myBestTC.modresp{1}{3,3}(:,1),'r')










            
            
            
            
