function explore_data()
% explore_data  - look through dat
%
% make an interface to look through data
%
% ds 2022-03-29

%% set up window

% obj = onCleanup(@()close('all'));

% make display window big:
screensize = get( groot(), 'Screensize' );

w.h = figure(99);
set(w.h, 'name', 'Display window for pRF data');
set(w.h, 'position', [10, 10, 0.75*screensize(3), 0.75*screensize(3)]);

MAXN = 128;

models = {'_2ddouble_', ...
    '_1d_', ...
    '_1dt_', ...
    '_sixteen_'};

cmaps = getColors('cmaps.txt');
cmaps = cmaps{1}; % make into correct shape

% this where we get the data
global data;
if isempty(data) % and various other checks // cache flag??
    grabData();
end
% provides data, mysubs, nModels, datafilenames, etc....


%% pass in info about display fig
localData.data = data;
localData.w = w;
localData.models = models; % for turning strings to iModel

paramsInfo = {...
    {'montage',0,'type=checkbox','callback', @displayCallback,'callbackArg', localData, 'passParams=1', 'Show in montage mode?'},...
    {'model',models,'callback', @displayCallback,'callbackArg', localData, 'passParams=1', 'Which of the 4 models to display'},... % will be resolved in callback
    {'useMean',0,'type=checkbox','callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'doPCA',1,'type=checkbox','callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'N',16,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %d]', MAXN),'callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'iSub',1,'incdec=[-1 1]','round=1','minmax=[1 10]','callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'iDigit',1,'incdec=[-1 1]','round=1','minmax=[1 4]','callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'iBaa',1,'incdec=[-1 1]','round=1','minmax=[1 4]','callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'colormap',cmaps,'callback', @displayCallback,'callbackArg', localData, 'passParams=1'},...
    {'colorbar',1,'type=checkbox','callback', @displayCallback,'callbackArg', localData, 'passParams=1'}...
    };

%%
p = mrParamsDialog(paramsInfo,'Controls',[]);

% close display windos / try/catch also possible
try
    close(w.h)
catch
    disp('window already gone...')
end
fprintf('BYE!\n')

end

function displayCallback(localData, params)

% unpack params that get passed in
iSub = params.iSub;
iDigit = params.iDigit;
iBaa = params.iBaa;

%iModel = strmatch(params.model,localData.models); % find index (tells you whether 1,2, 3 ...);
iModel = find(strcmpi(params.model,localData.models));

% factor out getting the data
% this function uses different model_* functions to return RF shapes
[thisrf, thisr2, d_desc] = return_rfs(localData.data, iSub, iDigit, iBaa, iModel, params.N);

% how about instead of an RF
% we try a PCA approach
thisrf_reshaped = reshape(thisrf, [size(thisrf,1).*size(thisrf,1) size(thisrf,3)]);
% need to flip it as ROWS = Observations, and COLUMNS = Variables;
% so e.g. we have 48 observations (voxels), and 16 variables (each pRF
% shape), for one ROI, e.g. D2, PD2
Xt = transpose(thisrf_reshaped); % - need to check this w/Denis
[coeff, score, latent, tsquared, explained] = pca(Xt);
% each column of coeff is a principal component in descending order of
% component variance
topPC = coeff(:,1);
topPC_rf = reshape(topPC,4,4); % this is your principal component shape from top N prfs


f_ = figure(localData.w.h); % use the display figure

title_txt = sprintf('digit: %d, pd %d, m: %s', d_desc.currentDigit, ...
    d_desc.currentBaa, d_desc.currentModel);


if params.montage
    sz = ceil(sqrt(params.N));
    t_ = tiledlayout(sz, sz);
    % make sure we don't ask for data that isn't there (no voxelsin ROI)
    maxT = min(size(thisrf, 3), params.N);
    for iT = 1:maxT
        nexttile()
        imagesc(thisrf(:,:,iT));
        h = plot_contours(thisrf(:,:,iT), [], 'color','w', 'linewidth',2);
        pbaspect([1 1 1]);
        axis('xy')
        camroll(-90); %rotate?
    end
    title(t_, title_txt, 'interpreter','none');
    % force size to keep window consistent
    % montage(permute(thisrf, [2 1 3]),  'BorderSize', 10, 'Size', [sz sz]);
elseif params.useMean
    imagesc(mean(thisrf,3));
    h = plot_contours(thisrf, [], 'color','w', 'linewidth',2);
    pbaspect([1 1 1]);
    axis('ij')
elseif params.doPCA
    imagesc(topPC_rf)
    pbaspect([1 1 1]);
    text(1,1,[num2str(explained(1)) ' % explained'],'Color','white','FontSize',14')
%     pca_txt = sprintf('digit: %d, pd %d, m: %s', d_desc.currentDigit, ...
%     d_desc.currentBaa, d_desc.currentModel);
    
    
else
    imagesc(thisrf(:,:,end));
    h = plot_contours(thisrf(:,:,end), [], 'color','w', 'linewidth',2);
    
    pbaspect([1 1 1]);
    axis('ij')
end
title(title_txt, 'interpreter','none');

colormapfunc = str2func(params.colormap);
colormap(f_, colormapfunc() )



disp('-- displayCallback() callback');

end

function myCallback(src, evt)

disp('-- myCallback() callback');

end

function C = getColors(fname)

fid = fopen(fname);
C = textscan(fid, '%s');
fclose(fid);

end

function h = plot_contours(d, theLevel, varargin)
% plot_contours - make a plot of all RF contours
% spread
%
% ds 2022-03

% plot peak(s)
plotCentres = true;

% define a color map
if ieNotDefined('theLevel')
    theLevel = 0.9; % RF contour at this height!
end

% varargin for linespec

% loop through all of the d's in the 3rd dim... and get contours

nPRFs = size(d,3);
hold('on')
allX = [];
allY = [];
h_=[];
s_ = [];
for iPRF = 1:nPRFs
    
    currentD = d(:,:,iPRF);
    minD = min(currentD(:));
    rangeD = range(currentD(:));
    halfD = minD + rangeD/2;
    if halfD > 0
        % separate comp and plotting
        [~, h_] = contour(currentD, [halfD halfD], varargin{:}); % see docs for contour
        
        % plot max of RF centers
        [x0, y0, v] = find2dMax(currentD);
        allX = [allX; x0];
        allY = [allY; y0];
        
    else
        fprintf('(%d/%d) halfD: %.2f, minD: %.2f\n', iPRF, nPRFs, halfD, minD);
    end
end


if plotCentres
    s_= scatter(allX,allY,'filled', 'markerfacecolor', 'w', ...
        'markeredgecolor', 'k')
end

% output arg
h = {h_,s_};

end

function [x0, y0, v] = find2dMax(M)
% find2dMax - find location and value of max in 2d
%

[v,I] = max(M(:));
[y0,x0] = ind2sub(size(M),I);

end


