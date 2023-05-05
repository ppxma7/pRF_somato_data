%% load data
% Compare R2 and AIC
%
% Michael Asghar March 2023
%
% See Also installAIC.m
clear  variables
close all
mythresh = 0;
doFac = 0;
mynbins = 4;
mydir = '/Volumes/nemosine/prfsomato_fits_aic_march2023/';
subjectList = {'pRF1','pRF2','pRF3','pRF4','pRF6','pRF7','pRF8','pRF9','pRF10','pRF11','pRF12'};
%subjectList = {'pRF7'};


fudge = struct;

for ii = 1:length(subjectList)

    load([mydir subjectList{ii} '/prf_1d_nodiag2_nov.mat'],'prf_overlays')
    %load([mydir subjectList{ii} '/AIC_no2p_prf_1d_nodiag2_nov.mat'],'prf_overlays')
    r2_1d = prf_overlays(:,1);
    aic_1d = prf_overlays(:,6);

    %load([mydir subjectList{ii} '/AIC_no2p_prf_1dt_nodiag2_nov.mat'],'prf_overlays')
    load([mydir subjectList{ii} '/prf_1dt_nodiag2_nov.mat'],'prf_overlays')
    r2_1dt = prf_overlays(:,1);
    aic_1dt = prf_overlays(:,6);

    %load([mydir subjectList{ii} '/AIC_no2p_prf_sixteen_nodiag2_com_nov.mat'],'prf_overlays')
    load([mydir subjectList{ii} '/prf_sixteen_nodiag2_com_nov.mat'],'prf_overlays')
    r2_16 = prf_overlays(:,1);
    aic_16 = prf_overlays(:,9);

    %load([mydir subjectList{ii} '/AIC_no2p_prf_2ddouble_dd_nov.mat'],'prf_overlays')
    load([mydir subjectList{ii} '/prf_2ddouble_dd_nov.mat'],'prf_overlays')
    r2_2d = prf_overlays(:,1);
    aic_2d = prf_overlays(:,7);

    fudge(ii).name = subjectList{ii};

    fudge(ii).r2_1d = r2_1d;
    fudge(ii).r2_1dt = r2_1dt;
    fudge(ii).r2_16 = r2_16;
    fudge(ii).r2_2d = r2_2d;

    fudge(ii).aic_1d = aic_1d;
    fudge(ii).aic_1dt = aic_1dt;
    fudge(ii).aic_16 = aic_16;
    fudge(ii).aic_2d = aic_2d;

end


% unpack

if length(subjectList)>1
    r2_1d        = [fudge(1).r2_1d; fudge(2).r2_1d; fudge(3).r2_1d; fudge(4).r2_1d; fudge(5).r2_1d; fudge(6).r2_1d; fudge(7).r2_1d; fudge(8).r2_1d; fudge(9).r2_1d; fudge(10).r2_1d; fudge(11).r2_1d  ];
    r2_1dt       = [fudge(1).r2_1dt;fudge(2).r2_1dt;fudge(3).r2_1dt;fudge(4).r2_1dt;fudge(5).r2_1dt;fudge(6).r2_1dt;fudge(7).r2_1dt;fudge(8).r2_1dt;fudge(9).r2_1dt;fudge(10).r2_1dt;fudge(11).r2_1dt ];
    r2_16        = [fudge(1).r2_16; fudge(2).r2_16; fudge(3).r2_16; fudge(4).r2_16; fudge(5).r2_16; fudge(6).r2_16; fudge(7).r2_16; fudge(8).r2_16; fudge(9).r2_16; fudge(10).r2_16; fudge(11).r2_16  ];
    r2_2d        = [fudge(1).r2_2d; fudge(2).r2_2d; fudge(3).r2_2d; fudge(4).r2_2d; fudge(5).r2_2d; fudge(6).r2_2d; fudge(7).r2_2d; fudge(8).r2_2d; fudge(9).r2_2d; fudge(10).r2_2d; fudge(11).r2_2d  ];

    aic_1d        = [fudge(1).aic_1d; fudge(2).aic_1d; fudge(3).aic_1d; fudge(4).aic_1d; fudge(5).aic_1d; fudge(6).aic_1d; fudge(7).aic_1d; fudge(8).aic_1d; fudge(9).aic_1d; fudge(10).aic_1d; fudge(11).aic_1d  ];
    aic_1dt       = [fudge(1).aic_1dt;fudge(2).aic_1dt;fudge(3).aic_1dt;fudge(4).aic_1dt;fudge(5).aic_1dt;fudge(6).aic_1dt;fudge(7).aic_1dt;fudge(8).aic_1dt;fudge(9).aic_1dt;fudge(10).aic_1dt;fudge(11).aic_1dt ];
    aic_16        = [fudge(1).aic_16; fudge(2).aic_16; fudge(3).aic_16; fudge(4).aic_16; fudge(5).aic_16; fudge(6).aic_16; fudge(7).aic_16; fudge(8).aic_16; fudge(9).aic_16; fudge(10).aic_16; fudge(11).aic_16  ];
    aic_2d        = [fudge(1).aic_2d; fudge(2).aic_2d; fudge(3).aic_2d; fudge(4).aic_2d; fudge(5).aic_2d; fudge(6).aic_2d; fudge(7).aic_2d; fudge(8).aic_2d; fudge(9).aic_2d; fudge(10).aic_2d; fudge(11).aic_2d  ];
else

    r2_1d = fudge(1).r2_1d;
    r2_1dt = fudge(1).r2_1dt;
    r2_16 = fudge(1).r2_16;
    r2_2d = fudge(1).r2_2d;
    aic_1d = fudge(1).aic_1d;
    aic_1dt = fudge(1).aic_1dt;
    aic_16 = fudge(1).aic_16;
    aic_2d = fudge(1).aic_2d;


end


% remove nans
myidx_1d = ~isnan(r2_1d);
myidx_1dt = ~isnan(r2_1dt);
myidx_16 = ~isnan(r2_16);
myidx_2d = ~isnan(r2_2d);

r2_1d = r2_1d(myidx_1d);
r2_1dt = r2_1dt(myidx_1dt);
r2_16 = r2_16(myidx_16);
r2_2d = r2_2d(myidx_2d);

aic_1d = aic_1d(myidx_1d);
aic_1dt = aic_1dt(myidx_1dt);
aic_16 = aic_16(myidx_16);
aic_2d = aic_2d(myidx_2d);

% remove zeros - but this leads to different models having different
% lengths!

% myidxzero_1d = r2_1d>0;
% myidxzero_1dt = r2_1dt>0;
% myidxzero_16 = r2_16>0;
% myidxzero_2d = r2_2d>0;
% 
% r2_1d = r2_1d(myidxzero_1d);
% r2_1dt = r2_1dt(myidxzero_1dt);
% r2_16 = r2_16(myidxzero_16);
% r2_2d = r2_2d(myidxzero_2d);
% 
% aic_1d = aic_1d(myidxzero_1d);
% aic_1dt = aic_1dt(myidxzero_1dt);
% aic_16 = aic_16(myidxzero_16);
% aic_2d = aic_2d(myidxzero_2d);

%%
% first plot mean R2

figure
allR2 = [r2_1d; r2_1dt; r2_16; r2_2d];
allAIC = [aic_1d; aic_1dt; aic_16; aic_2d];
mycols = [repmat({'1D Gaussian WD Model'},length(r2_1d),1);...
    repmat({'1D Gaussian BD Model'},length(r2_1d),1);...
    repmat({'Unconstrained Model'},length(r2_1d),1);...
    repmat({'2D Gaussian Digit Model'},length(r2_1d),1)];
thexx = repmat({'Models'},length(r2_1d),4);

g = gramm('x',thexx, 'y', allR2,'color',mycols);
g.stat_summary('geom',{'bar','black_errorbar'});
%g.geom_point()
g.set_names('x',[],'y','R2')
g.set_text_options('Font','Helvetica', 'base_size', 12)
g.set_point_options('base_size',12)
%g.axe_property('XGrid','on','YGrid','on','DataAspectRatio',[1 1 1])
g.set_order_options('x',0)
g.axe_property('XGrid','on','YGrid','on','YLim',[0 0.2])
%g.set_color_options('map',mycmap)

g.draw

filename1 = 'modelR2';
% g.export('file_name',filename1, ...
%     'export_path',...
%     '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
%     'file_type','pdf')




%% raincloud
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cl(1,:) = cb(4,:);
cl(2,:) = cb(5,:);
cl(3,:) = cb(6,:);
cl(4,:) = cb(7,:);
cl(5,:) = cb(8,:);
cl(6,:) = cb(9,:);

%myrow = 2;
%mycol = 2;
myalpha = 0.5;
myboxdodge = 0.1:0.2:3;
mydotdodge = 0.1:0.2:3;


figure %('Position', [100 100 1800 800])
h1 = raincloud_plot(r2_1d, 'box_on',1, 'color', cl(1,:),...
    'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(1), 'dot_dodge_amount', mydotdodge(1));
h2 = raincloud_plot(r2_1dt, 'box_on',1, 'color', cl(2,:),...
    'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(2), 'dot_dodge_amount', mydotdodge(2));
h3 = raincloud_plot(r2_16, 'box_on',1, 'color', cl(3,:),...
    'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(3), 'dot_dodge_amount', mydotdodge(3));
h4 = raincloud_plot(r2_2d, 'box_on',1, 'color', cl(4,:),...
    'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(4), 'dot_dodge_amount', mydotdodge(4));
% h4 = raincloud_plot(r24, 'box_on',1, 'color', cl(4,:),...
%     'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(4), 'dot_dodge_amount', mydotdodge(4));
% h5 = raincloud_plot(r25, 'box_on',1, 'color', cl(5,:),...
%     'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(5), 'dot_dodge_amount', mydotdodge(5));
% h6 = raincloud_plot(r26, 'box_on',1, 'color', cl(6,:),...
%     'alpha', myalpha, 'box_dodge',1, 'box_dodge_amount',myboxdodge(6), 'dot_dodge_amount', mydotdodge(6));
% legend([h1{1} h2{1} h3{1} h4{1} h5{1} h6{1}], {'HCTW', 'HCM', 'PostTW', 'PostM', 'PreTW', 'PreM'});
legend([h1{1} h2{1} h3{1} h4{1}],{'1D Gaussian WD','1D Gaussian BD','Unconstrained','2D Gaussian'});
ylim([-10 20])

set(gca,'FontSize',12)

filename2 = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/','raincloud2');
%print(filename2,'-dpdf')




%% now look at AIC

% matrix!
% what percentage of modelA vs modelB has lower AIC in modelA



prc(1) = sum(aic_1d<aic_1d)./length(aic_1d) .* 100;
prc(2) = sum(aic_1d<aic_1dt)./length(aic_1d) .* 100;
prc(3) = sum(aic_1d<aic_2d)./length(aic_1d) .* 100;
prc(4) = sum(aic_1d<aic_16)./length(aic_1d) .* 100;
prc(5) = sum(aic_1dt<aic_1d)./length(aic_1dt) .* 100;
prc(6) = sum(aic_1dt<aic_1dt)./length(aic_1dt) .* 100;
prc(7) = sum(aic_1dt<aic_2d)./length(aic_1dt) .* 100;
prc(8) = sum(aic_1dt<aic_16)./length(aic_1dt) .* 100;
prc(9) = sum(aic_2d<aic_1d)./length(aic_2d) .* 100;
prc(10) = sum(aic_2d<aic_1dt)./length(aic_2d) .* 100;
prc(11) = sum(aic_2d<aic_2d)./length(aic_2d) .* 100;
prc(12) = sum(aic_2d<aic_16)./length(aic_2d) .* 100;
prc(13) = sum(aic_16<aic_1d)./length(aic_16) .* 100;
prc(14) = sum(aic_16<aic_1dt)./length(aic_16) .* 100;
prc(15) = sum(aic_16<aic_2d)./length(aic_16) .* 100;
prc(16) = sum(aic_16<aic_16)./length(aic_16) .* 100;

mymat = transpose(reshape(prc,4,4));


figure, imagesc(mymat)
colormap viridis
colorbar
%clim([20 65])

xticks(1:4)
yticks(1:4)
xticklabels({'1D Gaussian WD','1D Gaussian BD','2D Gaussian','Unconstrained'});
yticklabels({'1D Gaussian WD','1D Gaussian BD','2D Gaussian','Unconstrained'});
xtickangle(45)
ytickangle(45)
ax = gca;
ax.FontSize = 12;
axis square
filename6 = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/','aicmatrix');
print(filename6,'-dpdf')
%% mean and std of raw AIC
% we want to do the prc diff calculation within subjects, then get a mean
% prc and std prc across subjects, so will get a prc X subject matrix

% clear nans first individually
toffee = fudge;
for tt = 1:length(subjectList)
    thedex_1d = ~isnan(toffee(tt).aic_1d);
    thedex_1dt = ~isnan(toffee(tt).aic_1dt);
    thedex_2d = ~isnan(toffee(tt).aic_2d);
    thedex_16 = ~isnan(toffee(tt).aic_16);
    toffee(tt).aic_1d = toffee(tt).aic_1d(thedex_1d);
    toffee(tt).aic_1dt = toffee(tt).aic_1dt(thedex_1dt);
    toffee(tt).aic_2d = toffee(tt).aic_2d(thedex_2d);
    toffee(tt).aic_16 = toffee(tt).aic_16(thedex_16);
end

for ff = 1:length(subjectList)
    prc_mstd(1,ff) = sum(toffee(ff).aic_1d<toffee(ff).aic_1d)./length(toffee(ff).aic_1d) .* 100;
    prc_mstd(2,ff) = sum(toffee(ff).aic_1d<toffee(ff).aic_1dt)./length(toffee(ff).aic_1d) .* 100;
    prc_mstd(3,ff) = sum(toffee(ff).aic_1d<toffee(ff).aic_2d)./length(toffee(ff).aic_1d) .* 100;
    prc_mstd(4,ff) = sum(toffee(ff).aic_1d<toffee(ff).aic_16)./length(toffee(ff).aic_1d) .* 100;
    prc_mstd(5,ff) = sum(toffee(ff).aic_1dt<toffee(ff).aic_1d)./length(toffee(ff).aic_1dt) .* 100;
    prc_mstd(6,ff) = sum(toffee(ff).aic_1dt<toffee(ff).aic_1dt)./length(toffee(ff).aic_1dt) .* 100;
    prc_mstd(7,ff) = sum(toffee(ff).aic_1dt<toffee(ff).aic_2d)./length(toffee(ff).aic_1dt) .* 100;
    prc_mstd(8,ff) = sum(toffee(ff).aic_1dt<toffee(ff).aic_16)./length(toffee(ff).aic_1dt) .* 100;
    prc_mstd(9,ff) = sum(toffee(ff).aic_2d<toffee(ff).aic_1d)./length(toffee(ff).aic_2d) .* 100;
    prc_mstd(10,ff) = sum(toffee(ff).aic_2d<toffee(ff).aic_1dt)./length(toffee(ff).aic_2d) .* 100;
    prc_mstd(11,ff) = sum(toffee(ff).aic_2d<toffee(ff).aic_2d)./length(toffee(ff).aic_2d) .* 100;
    prc_mstd(12,ff) = sum(toffee(ff).aic_2d<toffee(ff).aic_16)./length(toffee(ff).aic_2d) .* 100;
    prc_mstd(13,ff) = sum(toffee(ff).aic_16<toffee(ff).aic_1d)./length(toffee(ff).aic_16) .* 100;
    prc_mstd(14,ff) = sum(toffee(ff).aic_16<toffee(ff).aic_1dt)./length(toffee(ff).aic_16) .* 100;
    prc_mstd(15,ff) = sum(toffee(ff).aic_16<toffee(ff).aic_2d)./length(toffee(ff).aic_16) .* 100;
    prc_mstd(16,ff) = sum(toffee(ff).aic_16<toffee(ff).aic_16)./length(toffee(ff).aic_16) .* 100;
end

% now reshape to subjectList * 4x4
for ff = 1:length(subjectList)
    prc_mstd_mat(:,:,ff) = transpose(reshape(prc_mstd(:,ff),4,4));
end


mymeanmat = mean(prc_mstd_mat,3);
mystdmat = std(prc_mstd_mat,[],3);


figure('Position',[100 100 800 400])
tiledlayout(1,2)
nexttile
imagesc(mymeanmat)
colormap viridis
colorbar
clim([0 100])
xticks(1:4)
yticks(1:4)
xticklabels({'1D Gaussian WD','1D Gaussian BD','2D Gaussian','Unconstrained'});
yticklabels({'1D Gaussian WD','1D Gaussian BD','2D Gaussian','Unconstrained'});
xtickangle(45)
ytickangle(45)
ax = gca;
ax.FontSize = 12;
title('Mean')
axis square

nexttile
imagesc(mystdmat)
colormap viridis
colorbar
clim([0 100])
xticks(1:4)
yticks(1:4)
xticklabels({'1D Gaussian WD','1D Gaussian BD','2D Gaussian','Unconstrained'});
yticklabels({'1D Gaussian WD','1D Gaussian BD','2D Gaussian','Unconstrained'});
xtickangle(45)
ytickangle(45)
ax = gca;
ax.FontSize = 12;
title('Standard deviation')
axis square
 
filenamemeanstd = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/','aicmatrix_meanstd.png');
print(filenamemeanstd,'-dpng')

%% this is testing if there is a sig diff of AIC between models, across all subs
% concat all subs' AIC into 4 columns, one for each model
bla1d = {toffee.aic_1d}';
blav1d = vertcat(bla1d{:});

bla1dt = {toffee.aic_1dt}';
blav1dt = vertcat(bla1dt{:});

bla16 = {toffee.aic_16}';
blav16 = vertcat(bla16{:});

bla2d = {toffee.aic_2d}';
blav2d = vertcat(bla2d{:});

blavstack = [blav1d, blav1dt, blav16, blav2d];
[P,ANOVATAB,STATS] = anova1(blavstack);
[COMPARISON,MEANS,H,GNAMES] = multcompare(STATS);

%% this is testing if there is a sig diff of percentage differences, between subjects
% prc_mstd is the % diffs where each column is a sub

[P,ANOVATAB,STATS] = anova1(prc_mstd);
[COMPARISON,MEANS,H,GNAMES] = multcompare(STATS);


%% this is testing if there is a sig diff between the percentage differences

prc_mstd_flipped = transpose(prc_mstd);
prc_mstd_flipped(:,1) = [];
prc_mstd_flipped(:,5) = [];
prc_mstd_flipped(:,9) = [];
prc_mstd_flipped(:,13) = [];
[P,ANOVATAB,STATS] = anova1(prc_mstd_flipped);
[COMPARISON,MEANS,H,GNAMES] = multcompare(STATS);

compdex = COMPARISON(:,6)<0.05;
sigcomps = COMPARISON(compdex,:);

tbl = array2table(sigcomps,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=GNAMES(tbl.("Group A"));
tbl.("Group B")=GNAMES(tbl.("Group B"));

filenamemeanH = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/','aic_anova1_tbl');
filenamemeanHH = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/','aic_multc_tbl');
writecell(ANOVATAB,filenamemeanH,'FileType','spreadsheet')
writetable(tbl,filenamemeanHH,'FileType','spreadsheet')





%%

%
%figure('Position',[100 100 600 1200])
figure
tiledlayout(4,4)
col1 = [0.8500 0.3250 0.0980]; % orange
col2 = [0 0.4470 0.7410]; % blue
edgealpha = 0.2;
facealpha = 0.4;


nexttile
set(gca,'XTick',[], 'YTick', [])


nexttile
histogram(aic_1d,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_1dt,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_1d,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_2d,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_1d,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_16,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_1dt,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_1d,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
set(gca,'XTick',[], 'YTick', [])


nexttile
histogram(aic_1dt,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_2d,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_1dt,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_16,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_2d,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_1d,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_2d,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_1dt,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile

set(gca,'XTick',[], 'YTick', [])

nexttile
histogram(aic_2d,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_16,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_16,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_1d,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_16,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_1dt,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile
histogram(aic_16,'FaceColor',col1,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col1)
hold on
histogram(aic_2d,'FaceColor',col2,'EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'EdgeColor',col2)
axis tight

nexttile

set(gca,'XTick',[], 'YTick', [])

filename4 = fullfile('/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/','aichistmatrix_no2p.pdf');
print(filename4,'-dpdf')












