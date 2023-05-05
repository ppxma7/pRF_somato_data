%% try to refresh prfstat_testingBaseTips % prfstat_testingBAll
% April 2022
% ma
%
close all
clear variables
close all hidden

% here the beginning is the same as prfsummary2.m
thisorange = [0.85 0.325 0.098];
thisyellow = [0.929 0.694 0.125];
thisgreen = [0.466 0.6740 0.188];
thisblue = [0 0.447 0.741];

thiswhite = [1 1 1];
mycmap = [thisorange; thisyellow; thisgreen; thisblue];


mypath='/Volumes/nemosine/prfsomato_fits';
%mypath = '/Users/lpzds1/Library/CloudStorage/OneDrive-SharedLibraries-TheUniversityofNottingham/Michael_Sue - pRF_paper_draft/prfsomato_fits';

% range of data
% SUBJECTS
mysubs = {'prf1','prf2','prf3','prf4','prf6','prf7','prf8','prf9','prf10','prf11'};
% test case
% mysubs = {'prf1','prf2'};
nsubs = 10;
%howmanyrows = 4;

thisModel = 2;
doBa = 1;
doBasetips = 0;
aa = 3.5;
if doBasetips == 1
    howmanyrows = 4;
elseif doBa == 1
    howmanyrows = 4;
end

% DIGITS / PD locations
digits = 2:5;
baa = 2:howmanyrows+1;
mynbins = 4;
mycutoff = 48;
upperthresh = 15;

datafilenames = {'prf_2ddouble_dd_nov', ...
    'prf_1d_nodiag2_nov', ...
    'prf_1dt_nodiag2_nov', ...
    'prf_sixteen_nodiag2_com_nov'};

nModels = numel(datafilenames);

data = {}; % empty struct
results = {};

for iSub = 1:nsubs%:length(mysubs)
    
    currentSub = mysubs{iSub};
    
    tic()
    % for each data file
    for iModel = 1:nModels
        currentModel = datafilenames{iModel};
        
        % Read in all BA ROIs and indices
        for iBaa = 1:howmanyrows
            currentBaa = baa(iBaa);
            
            for iDigit = 1:4
                %currentBaa = baa(iBaa);
                currentDigit = digits(iDigit);
                % just pick one model (ROI indeces are all the same
                % load([mypath mysubs{iSub} '/baa/' 'D' num2str(iDig) 'pd' num2str(iBaa) '_16.mat'],'prf_overlays');
                currentROIFilename = [mypath filesep() currentSub filesep() '/baa/' sprintf('D%dpd%d_16.mat', currentDigit, currentBaa)];
                data{iSub}.roiIndex(iBaa, iDigit) = load(currentROIFilename, 'idx');
            end
        end
        
        %         % find the limits of the big roi
        %         % this is important, becuase it uses the BA limits to contain the
        %         % PD rois, otherwise we are involving responses from BA4, and 6
        %         % etc. because the big PD ROI extent is bigger than all 16 BA areas
        %         % combined
        %         A = struct2cell(data{iSub}.roiIndex);
        %         A = A(:);
        %         Astack = cell2mat(A);
        
        % basetips
        currentFilename = [mypath filesep() currentSub filesep() 'basetips' filesep() currentModel ];
        
        data{iSub}.model(iModel).data = load(currentFilename);
        data{iSub}.model(iModel).modelname = currentModel;
        
        % common, big, small
        %         [commonVoxels_up, iInData_up, iInROI_up] = intersect(data{iSub}.model(iModel).data.idx, Astack);
        %         data{iSub}.model(iModel).data.prf_overlays_limited = data{iSub}.model(iModel).data.prf_overlays(iInData_up,:);
        %         data{iSub}.model(iModel).data.idx_limited = commonVoxels_up;
        %
        data{iSub}.model(iModel).data.prf_overlays_limited = data{iSub}.model(iModel).data.prf_overlays;
        data{iSub}.model(iModel).data.idx_limited = data{iSub}.model(iModel).data.idx;
        
        % remake the converted params for PD
        thisDig = data{iSub}.model(iModel).data.prf_overlays_limited(:,2);
        thisPD = data{iSub}.model(iModel).data.prf_overlays_limited(:,3);
        
        thisDig_conv = convert2discretebins(thisDig,mynbins);
        thisPD_conv = convert2discretebins(thisPD,mynbins);
        
        data{iSub}.model(iModel).convDig = thisDig_conv;
        %data{iSub}.model(iModel).convPD = thisPD_conv;
        
        if iModel == 1 %2DGaussian
            
            % fudge factor as maps are flipped from pRFsomatoFit output for
            % PD
            pd_2d_1dex = thisPD_conv==1;
            pd_2d_2dex = thisPD_conv==2;
            pd_2d_3dex = thisPD_conv==3;
            pd_2d_4dex = thisPD_conv==4;
            thisPD_conv(pd_2d_1dex) = 4;
            thisPD_conv(pd_2d_2dex) = 3;
            thisPD_conv(pd_2d_3dex) = 2;
            thisPD_conv(pd_2d_4dex) = 1;
            
            data{iSub}.model(iModel).convPD = thisPD_conv;
        elseif iModel == 2 || iModel == 3
            % there is an abs situation here in pRFsomatoFit see line
        
            pd_2d_1dex = thisPD_conv==1;
            pd_2d_2dex = thisPD_conv==2;
            pd_2d_3dex = thisPD_conv==3;
            pd_2d_4dex = thisPD_conv==4;
            thisPD_conv(pd_2d_1dex) = (thisPD_conv(pd_2d_1dex).*-1)+5;
            thisPD_conv(pd_2d_2dex) = (thisPD_conv(pd_2d_2dex).*-1)+5;
            thisPD_conv(pd_2d_3dex) = (thisPD_conv(pd_2d_3dex).*-1)+5;
            thisPD_conv(pd_2d_4dex) = (thisPD_conv(pd_2d_4dex).*-1)+5;
            
            data{iSub}.model(iModel).convPD = thisPD_conv;
            
            
        else
            data{iSub}.model(iModel).convPD = thisPD_conv;
        end
        
        % now we can overlap and index to get PD ROIs.
        for iPD = 1:4
            for iDigit = 1:4
                thisDex = data{iSub}.model(iModel).convDig == iDigit & data{iSub}.model(iModel).convPD == iPD;
                thisROI = data{iSub}.model(iModel).data.idx(thisDex);
                % fudge to fit a struct into a cell array
                tmpStruct = struct();
                tmpStruct.idx = thisROI;
                data{iSub}.model(iModel).roiIndexPD(iPD,iDigit) = tmpStruct;
                
            end
        end
        
    end
    
    
    
    
    t = toc();
    fprintf('took: %.2fs to read sub %d\n', t, iSub);
    
    
end

%%
% so at this point we have 2 sets of ROI indices.
% 1 4x4 struct for the BAA based on the ROI loading
% 1 4x4 struct for the PD based on how I did some mixing on PD values
% above.


% look at the 'ex' field for this information
if thisModel == 1
    myadjr = 6;
    rfx = 4;
    rfy = 5;
elseif thisModel == 2 || thisModel == 3
    myadjr = 5;
    rfx = 4;
elseif thisModel == 4
    myadjr = 8;
    rfx = 7;
end



if doBa == 1
    sep = 'BAA';
elseif doBasetips == 1
    sep = 'BASETIPS';
end

for iSub = 1:nsubs%:10
    for iLoc = 1:howmanyrows
        for iDigit = 1:4
            
            currentSub = mysubs{iSub};
            currentDigit = digits(iDigit);
            currentBaa = baa(iLoc);
            currentModel = datafilenames{thisModel};
            currentIDX = data{iSub}.model(iModel).data.idx_limited;
            % current data
            currentD = data{iSub}.model(thisModel).data;
            
            currentO = data{iSub}.model(thisModel).data.prf_overlays_limited;
            
            % ROI - subset that's in the particular ROI -- same as in linearCoords
            % here is where we edit this for either BAA or PD
            if doBa == 1
                %currentROI = data{iSub}.roiIndex(iDigit, iBaa).idx;
                %currentO = data{iSub}.model(iModel).data.prf_overlays;
                currentROI = data{iSub}.roiIndex(iLoc, iDigit).idx;
                
            elseif doBasetips == 1
                %keyboard
                % need to create the idx here
                currentROI = data{iSub}.model(thisModel).roiIndexPD(iLoc, iDigit).idx;
                
            end
            
            theR2 = currentO(:,myadjr);
            
            % grab volumes
            BA_vols_no_model(iLoc,iDigit, iSub) = length(currentROI);
            
            % all DATA / params
            %extractedParams = currentD.params(:, ~invalidCoords);
            
            [commonVoxels, iInData, iInROI] = intersect(currentIDX, currentROI );
            
            %grab volumes of model specific
            BA_vols_models(iLoc,iDigit,iSub) = length(commonVoxels);
            
            roiR2 = theR2(iInData);
            roiR2 = roiR2(:);
            
            [topN, topNidx] = maxk(roiR2, mycutoff); % to be consistent with PCA later
            
            matchingLinearIdx = iInData(topNidx);
            %matchingLinearIdx = iInData;
            % okay we don't care about pRFs here, only just getting the
            % indices for each ROI, then we can index into pRF size values
            % and plot these.
            
            theRFx = currentO(:,rfx);
            roiRFx = theRFx(matchingLinearIdx);
            roiRFx = roiRFx(:);
            roiRFx(roiRFx>upperthresh)=[]; % clean up super high values
            if thisModel == 1
                theRFy = currentO(:,rfy);
                roiRFy = theRFy(matchingLinearIdx);
                roiRFy = roiRFy(:);
                roiRFy(roiRFy>upperthresh)=[]; % clean up super high values
                tmpStruct = struct();
                tmpStruct.rfy = roiRFy;
                results{iSub}.rfy(iLoc,iDigit) = tmpStruct;
            end
            
            tmpStruct = struct();
            tmpStruct.rfx = roiRFx;
            results{iSub}.rfx(iLoc,iDigit) = tmpStruct;
            
        end
    end 
end

%save(sprintf('BA_vols_%s_%s.mat',sep,currentModel),'BA_vols_models')

% double up for 1D
if thisModel == 2 || thisModel == 3
    thisModel = 3;
    for iSub = 1:nsubs
        for iLoc = 1:howmanyrows
            for iDigit = 1:4
                currentSub = mysubs{iSub};
                currentDigit = digits(iDigit);
                currentBaa = baa(iLoc);
                currentModel = datafilenames{thisModel};
                currentIDX = data{iSub}.model(iModel).data.idx_limited;
                currentD = data{iSub}.model(thisModel).data;
                currentO = data{iSub}.model(thisModel).data.prf_overlays_limited;
                if doBa == 1
                    currentROI = data{iSub}.roiIndex(iLoc, iDigit).idx;
                elseif doBasetips == 1
                    %currentROI = data{iSub}.roiIndexPD(iLoc, iDigit).idx;
                    currentROI = data{iSub}.model(thisModel).roiIndexPD(iLoc, iDigit).idx;
                end
                
                theR2 = currentO(:,myadjr);
                BA_vols_no_model2(iLoc,iDigit, iSub) = length(currentROI);
                
                [commonVoxels, iInData, iInROI] = intersect(currentIDX, currentROI );
                %grab volumes of ones that have activation
                BA_vols_models2(iLoc,iDigit,iSub) = length(commonVoxels);
                
                roiR2 = theR2(iInData);
                roiR2 = roiR2(:);
                
                [topN, topNidx] = maxk(roiR2, mycutoff); % to be consistent with PCA later
                
                matchingLinearIdx = iInData(topNidx);
                
                theRFx = currentO(:,rfx);
                roiRFx = theRFx(matchingLinearIdx);
                roiRFx = roiRFx(:);
                roiRFx(roiRFx>upperthresh)=[]; % clean up super high values
                
                tmpStruct = struct();
                tmpStruct.rfx = roiRFx;
                resultsT{iSub}.rfx(iLoc,iDigit) = tmpStruct;
                
            end
        end
    end
    
    
    
    
    %save(sprintf('BA_vols_%s_%s.mat',sep,currentModel),'BA_vols_models2')
    
    
    
    thisModel = 2;
    
    

end

%% now plot
% from here, it goes
% Column is a subject
% Rows: D2 D2 D2 D2 D3 etc.
%     :pd2 pd3 pd4 pd4 pd2 etc...

% unpack and average across subjects
for iSub = 1:nsubs
    tmp = struct2cell(results{iSub}.rfx);
    tmp = tmp(:);
    rfx_sizes{:,iSub} = tmp;
    %rfx_sizes(:,iSub) = rfx_sizes(:);
    rfx_means(:,iSub) = cellfun(@nanmean,tmp);
    
    if thisModel == 1
        tmpy = struct2cell(results{iSub}.rfy);
        tmpy = tmpy(:);
        rfy_sizes{:,iSub} = tmpy;
        rfy_means(:,iSub) = cellfun(@nanmean,tmpy);
    elseif thisModel == 2 || thisModel == 3
        rfy_sizes = struct2cell(resultsT{iSub}.rfx);
        rfy_sizes = rfy_sizes(:);
        rfy_means(:,iSub) = cellfun(@nanmean,rfy_sizes);
    end
end


% rfx_means = normalize(rfx_means,2,'scale');
% if thisModel ~=4
%     rfy_means = normalize(rfy_means,2,'scale');
% end


d2pd2 = rfx_means(1,:);
d2pd3 = rfx_means(2,:);
d2pd4 = rfx_means(3,:);
d2pd5 = rfx_means(4,:);

d3pd2 = rfx_means(5,:);
d3pd3 = rfx_means(6,:);
d3pd4 = rfx_means(7,:);
d3pd5 = rfx_means(8,:);

d4pd2 = rfx_means(9,:);
d4pd3 = rfx_means(10,:);
d4pd4 = rfx_means(11,:);
d4pd5 = rfx_means(12,:);

d5pd2 = rfx_means(13,:);
d5pd3 = rfx_means(14,:);
d5pd4 = rfx_means(15,:);
d5pd5 = rfx_means(16,:);

PD2m = [mean(d2pd2(:)); mean(d3pd2(:)); mean(d4pd2(:)); mean(d5pd2(:))];
PD3m = [mean(d2pd3(:)); mean(d3pd3(:)); mean(d4pd3(:)); mean(d5pd3(:))];
PD4m = [mean(d2pd4(:)); mean(d3pd4(:)); mean(d4pd4(:)); mean(d5pd4(:))];
PD5m = [mean(d2pd5(:)); mean(d3pd5(:)); mean(d4pd5(:)); mean(d5pd5(:))];

PD2mstd = [std(d2pd2(:))./sqrt(length(d2pd2)); std(d3pd2(:))./sqrt(length(d2pd2)); std(d4pd2(:))./sqrt(length(d2pd2)); std(d5pd2(:))./sqrt(length(d2pd2))];
PD3mstd = [std(d2pd3(:))./sqrt(length(d2pd2)); std(d3pd3(:))./sqrt(length(d2pd2)); std(d4pd3(:))./sqrt(length(d2pd2)); std(d5pd3(:))./sqrt(length(d2pd2))];
PD4mstd = [std(d2pd4(:))./sqrt(length(d2pd2)); std(d3pd4(:))./sqrt(length(d2pd2)); std(d4pd4(:))./sqrt(length(d2pd2)); std(d5pd4(:))./sqrt(length(d2pd2))];
PD5mstd = [std(d2pd5(:))./sqrt(length(d2pd2)); std(d3pd5(:))./sqrt(length(d2pd2)); std(d4pd5(:))./sqrt(length(d2pd2)); std(d5pd5(:))./sqrt(length(d2pd2))];

PD2 = [d2pd2(:); d3pd2(:); d4pd2(:); d5pd2(:)];
PD3 = [d2pd3(:); d3pd3(:); d4pd3(:); d5pd3(:)];
PD4 = [d2pd4(:); d3pd4(:); d4pd4(:); d5pd4(:)];
PD5 = [d2pd5(:); d3pd5(:); d4pd5(:); d5pd5(:)];

D2 = [d2pd2(:); d2pd3(:); d2pd4(:); d2pd5(:)];
D3 = [d3pd2(:); d3pd3(:); d3pd4(:); d3pd5(:)];
D4 = [d4pd2(:); d4pd3(:); d4pd4(:); d4pd5(:)];
D5 = [d5pd2(:); d5pd3(:); d5pd4(:); d5pd5(:)];

if thisModel == 1 || thisModel == 2 || thisModel == 3
    d2pd2y = rfy_means(1,:);
    d2pd3y = rfy_means(2,:);
    d2pd4y = rfy_means(3,:);
    d2pd5y = rfy_means(4,:);
    
    d3pd2y = rfy_means(5,:);
    d3pd3y = rfy_means(6,:);
    d3pd4y = rfy_means(7,:);
    d3pd5y = rfy_means(8,:);
    
    d4pd2y = rfy_means(9,:);
    d4pd3y = rfy_means(10,:);
    d4pd4y = rfy_means(11,:);
    d4pd5y = rfy_means(12,:);
    
    d5pd2y = rfy_means(13,:);
    d5pd3y = rfy_means(14,:);
    d5pd4y = rfy_means(15,:);
    d5pd5y = rfy_means(16,:);
    
    PD2ym = [mean(d2pd2y(:)); mean(d3pd2y(:)); mean(d4pd2y(:)); mean(d5pd2y(:))];
    PD3ym = [mean(d2pd3y(:)); mean(d3pd3y(:)); mean(d4pd3y(:)); mean(d5pd3y(:))];
    PD4ym = [mean(d2pd4y(:)); mean(d3pd4y(:)); mean(d4pd4y(:)); mean(d5pd4y(:))];
    PD5ym = [mean(d2pd5y(:)); mean(d3pd5y(:)); mean(d4pd5y(:)); mean(d5pd5y(:))];
    
    PD2ymstd = [std(d2pd2y(:))./sqrt(length(d2pd2)); std(d3pd2y(:))./sqrt(length(d2pd2)); std(d4pd2y(:))./sqrt(length(d2pd2)); std(d5pd2y(:))./sqrt(length(d2pd2))];
    PD3ymstd = [std(d2pd3y(:))./sqrt(length(d2pd2)); std(d3pd3y(:))./sqrt(length(d2pd2)); std(d4pd3y(:))./sqrt(length(d2pd2)); std(d5pd3y(:))./sqrt(length(d2pd2))];
    PD4ymstd = [std(d2pd4y(:))./sqrt(length(d2pd2)); std(d3pd4y(:))./sqrt(length(d2pd2)); std(d4pd4y(:))./sqrt(length(d2pd2)); std(d5pd4y(:))./sqrt(length(d2pd2))];
    PD5ymstd = [std(d2pd5y(:))./sqrt(length(d2pd2)); std(d3pd5y(:))./sqrt(length(d2pd2)); std(d4pd5y(:))./sqrt(length(d2pd2)); std(d5pd5y(:))./sqrt(length(d2pd2))];
    
    PD2y = [d2pd2y(:); d3pd2y(:); d4pd2y(:); d5pd2y(:)];
    PD3y = [d2pd3y(:); d3pd3y(:); d4pd3y(:); d5pd3y(:)];
    PD4y = [d2pd4y(:); d3pd4y(:); d4pd4y(:); d5pd4y(:)];
    PD5y = [d2pd5y(:); d3pd5y(:); d4pd5y(:); d5pd5y(:)];
    
    D2y = [d2pd2y(:); d2pd3y(:); d2pd4y(:); d2pd5y(:)];
    D3y = [d3pd2y(:); d3pd3y(:); d3pd4y(:); d3pd5y(:)];
    D4y = [d4pd2y(:); d4pd3y(:); d4pd4y(:); d4pd5y(:)];
    D5y = [d5pd2y(:); d5pd3y(:); d5pd4y(:); d5pd5y(:)];


end





%
DIGDEX = [repmat({'D2'},nsubs,1); repmat({'D3'},nsubs,1); repmat({'D4'},nsubs,1);repmat({'D5'},nsubs,1)];
DIGDEX2 = [repmat({'PD1'},nsubs.*4,1); repmat({'PD2'},nsubs.*4,1); repmat({'PD3'},nsubs.*4,1);repmat({'PD4'},nsubs.*4,1)];
DIGDEX3 = [repmat({'BA3a'},nsubs.*4,1); repmat({'BA3b'},nsubs.*4,1); repmat({'BA1'},nsubs.*4,1);repmat({'BA2'},nsubs.*4,1)];


%% explode plot
if thisModel == 1 || thisModel == 2 || thisModel == 3
    clear g
    figure('Position',[100 100 1400 400])
    g(1,1) = gramm('x',PD2y,'y',PD2,'color',DIGDEX);
    g(1,1).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,1).set_title('PD1')
    else
        g(1,1).set_title('BA3a')
    end
    g(1,2) = gramm('x',PD3y,'y',PD3,'color',DIGDEX);
    g(1,2).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,2).set_title('PD2')
    else
        g(1,2).set_title('BA3b')
    end
    g(1,3) = gramm('x',PD4y,'y',PD4,'color',DIGDEX);
    g(1,3).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,3).set_title('PD3')
    else
        g(1,3).set_title('BA1')
    end
    g(1,4) = gramm('x',PD5y,'y',PD5,'color',DIGDEX);
    g(1,4).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,4).set_title('PD4')
    else
        g(1,4).set_title('BA2')
    end
    
    
    
    g.set_names('x','Between Digit Direction','y','Within Digit Direction')
    g.set_text_options('Font','Helvetica', 'base_size', 16)
    g.set_point_options('base_size',10)
    g.axe_property('XGrid','on','YGrid','on','YLim',[0 5], 'XLim',[0 5],'DataAspectRatio',[1 1 1])
    g.set_order_options('x',0)
    g.set_color_options('map',mycmap)
    g.draw
    
    x = 0:10;
    y = 0:10;
    
    g(1,1).update('x',x,'y',y)
    g(1,1).geom_line
    g(1,2).update('x',x,'y',y)
    g(1,2).geom_line
    g(1,3).update('x',x,'y',y)
    g(1,3).geom_line
    g(1,4).update('x',x,'y',y)
    g(1,4).geom_line
    
    
    
    g.set_color_options('map',[0.5 0.5 0.5])
    g.draw()
    filename = sprintf(['michael_explode_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
    
    
else
    %DIGDEX3 = [{'PD2'},{'PD3'},{'PD4'},{'PD5'}];
    %DIGDEX4 = [repmat(DIGDEX3,4,1)];
    
    clear g
    figure('Position',[100 100 1400 400])
    g(1,1) = gramm('x',DIGDEX,'y',PD2,'color',DIGDEX);
    g(1,1).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,1).set_title('PD1')
    else
        g(1,1).set_title('BA3a')
    end
    g(1,2) = gramm('x',DIGDEX,'y',PD3,'color',DIGDEX);
    g(1,2).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,2).set_title('PD2')
    else
        g(1,2).set_title('BA3b')
    end
    g(1,3) = gramm('x',DIGDEX,'y',PD4,'color',DIGDEX);
    g(1,3).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,3).set_title('PD3')
    else
        g(1,3).set_title('BA1')
    end
    g(1,4) = gramm('x',DIGDEX,'y',PD5,'color',DIGDEX);
    g(1,4).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,4).set_title('PD4')
    else
        g(1,4).set_title('BA2')
    end
    
    
    
    
    g.set_names('x','Digit','y','Moment of Inertia')
    g.set_text_options('Font','Helvetica', 'base_size', 16)
    g.set_point_options('base_size',10)
    g.axe_property('XGrid','on','YGrid','on','YLim',[5 8],'DataAspectRatio',[1 1 1])
    g.set_order_options('x',0)
    g.set_color_options('map',mycmap)
    g.draw
    filename = sprintf(['michael_explode_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
    
    
end

%% now collapse above plot
 %axis limits
% stat_summary('geom',{'point' 'errorbar'},'dodge',0.3,'width',0.5);
DIGDEXm = [repmat({'D2'},1,1); repmat({'D3'},1,1); repmat({'D4'},1,1);repmat({'D5'},1,1)];
if strcmpi(sep,'BASETIPS')
    DIGDEX2m = [repmat({'PD1'},4,1); repmat({'PD2'},4,1); repmat({'PD3'},4,1);repmat({'PD4'},4,1)];
else
    DIGDEX2m = [repmat({'BA3a'},4,1); repmat({'BA3b'},4,1); repmat({'BA1'},4,1);repmat({'BA2'},4,1)];
    
    
end

mycolor1 = [thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue];



mycolorFace = [thisorange;thisorange;thiswhite;thiswhite;...
    thisyellow;thisyellow;thiswhite;thiswhite;...
    thisgreen;thisgreen;thiswhite;thiswhite;...
    thisblue;thisblue;thiswhite;thiswhite];

mymarker = ['s';'s';'s';'s';...
    '^';'^';'^';'^';...
    'o';'o';'o';'o';...
    'p';'p';'p';'p'];




if thisModel == 1 || thisModel == 2 || thisModel == 3
    clear g
    %figure('Position',[100 100 1400 400])
    figure
    
    g = gramm('x',[PD2ym;PD3ym;PD4ym;PD5ym],'y',[PD2m;PD3m;PD4m;PD5m],'color',[DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm],'marker',DIGDEX2m);
    
    %g.stat_summary('geom',{'point' 'errorbar'});
    g.geom_point()
    
    g.set_names('x','Between Digit Direction','y','Within Digit Direction')
    g.set_text_options('Font','Helvetica', 'base_size', 16)
    g.set_point_options('base_size',12)
    g.axe_property('XGrid','on','YGrid','on','YLim',[0 aa], 'XLim',[0 aa],'DataAspectRatio',[1 1 1])
    g.set_order_options('x',0)
    g.set_color_options('map',mycmap)
    g.draw
    
    x = 0:10;
    y = 0:10;
    
    g.update('x',x,'y',y)
    g.geom_line
    g.set_color_options('map',[0.5 0.5 0.5])
    g.draw
    
    filename = sprintf(['wang_explode_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
    % wang plot
    
    DYstack = [PD2ym;PD3ym;PD4ym;PD5ym];
    DXstack = [PD2m;PD3m;PD4m;PD5m];
    DYstack_std = [PD2ymstd;PD3ymstd;PD4ymstd;PD5ymstd];
    DXstack_std = [PD2mstd;PD3mstd;PD4mstd;PD5mstd];
    
    
    figure
    
    for kk = 1:16
        e = errorbar(DYstack(kk),DXstack(kk), DYstack_std(kk), DXstack_std(kk),...
            'both','o','MarkerEdgeColor', mycolor1(kk,:),'MarkerSize',10,...
            'MarkerFaceColor',mycolor1(kk,:),'Marker',mymarker(kk));
        e.Color = mycolor1(kk,:);
        hold on
    end
    
    axis square
    grid on
    x = 0:5;
    y = 0:5;
    plot(x,y,'k--')
    xlim([0 aa])
    ylim([0 aa])
    hold off
    xlabel('Between Digit Direction')
    ylabel('Within Digit Direction')
    %legend
    filename = sprintf(['/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/wang_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    print('-dpdf', filename)
    
    [H,P,CI,STATS] = ttest2(DYstack, DXstack);
    
elseif thisModel == 4
    clear g
    %figure('Position',[100 100 1400 400])
    figure
    
    g = gramm('x',DIGDEX2m,'y',[PD2m;PD3m;PD4m;PD5m],'color',[DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm]);
    
    %g.stat_summary('geom',{'point' 'errorbar'});
    g.geom_jitter()
    
    g.set_names('x','Digit','y','Moment of Inertia')
    g.set_text_options('Font','Helvetica', 'base_size', 16)
    g.set_point_options('base_size',12)
    g.axe_property('XGrid','on','YGrid','on','YLim',[5 8],'DataAspectRatio',[1 1 1])
    g.set_order_options('x',0)
    g.set_color_options('map',mycmap)
    g.draw
    
    filename = sprintf(['wang_explode_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
    
    
    DXstack = [PD2m;PD3m;PD4m;PD5m];
    DXstack_std = [PD2mstd;PD3mstd;PD4mstd;PD5mstd];
    bigcolor2 = [1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5,3.2,3.3,3.4,3.5,4.2,4.3,4.4,4.5];
    bigcolor3 = [1.35 2.35 3.35 4.35];
    
    
    figure()
    hold on
    
    for kk = 1:16
        e = errorbar(bigcolor2(kk),DXstack(kk), DXstack_std(kk),...
            'o','MarkerEdgeColor', mycolor1(kk,:),'MarkerSize',10,...
            'MarkerFaceColor',mycolor1(kk,:));
        e.Color = mycolor1(kk,:);
        hold on
    end
    xticks([1.35 2.35 3.35 4.35])
    xticklabels({'PD1','PD2','PD3','PD4'})
    xlim([1 4.5])
    
    
    
    
    
    
    %     %xticks([1.35 2.35 3.35 4.35])
    %     if strcmpi(sep,'BASETIPS')
    %        xticklabels({'PD2','PD3','PD4','PD5'})
    %     else
    %        %xticklabels({'PD2','PD3','PD4','PD5'})
    %        xticklabels({'BA3a','BA3b','BA1','BA2'})
    %     end
    %
    %xticklabels({'PD2','PD3','PD4','PD5'})
    %yticks([5:5:8])
    axis square
    grid on
    %xlim([1 4.5])
    ylim([5 8])
    hold off
    ylabel('Moment of Inertia')
    %xlabel('PD')
    filename = sprintf(['/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/wang_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    print('-dpdf', filename)
    
end


%% try all subjects wang plot
%
if thisModel ~= 4
clear g
figure
g = gramm('x',[PD2y;PD3y;PD4y;PD5y],'y',[PD2;PD3;PD4;PD5],'color',[DIGDEX;DIGDEX;DIGDEX;DIGDEX],'marker',DIGDEX2);
g.stat_summary('geom',{'point' 'errorbar'});

g.set_names('x','Between Digit Direction','y','Within Digit Direction')
g.set_text_options('Font','Helvetica', 'base_size', 16)
g.set_point_options('base_size',12)
g.axe_property('XGrid','on','YGrid','on','YLim',[0 5], 'XLim',[0 5],'DataAspectRatio',[1 1 1])
g.set_order_options('x',0)
g.set_color_options('map',mycmap)
g.draw

x = 0:10;
y = 0:10;

g.update('x',x,'y',y)
g.geom_line
g.set_color_options('map',[0.5 0.5 0.5])
g.draw

end

%% now sue plot

if strcmpi(sep,'BAA')
    whichSeparation = DIGDEX3; %BA
elseif strcmpi(sep,'BASETIPS')
    whichSeparation = DIGDEX2; %Basetips
end

if thisModel == 1 || thisModel == 2
    clear g
    figure('Position',[100 100 1000 500])
    g(1,1) = gramm('x',whichSeparation,'y',[PD2;PD3;PD4;PD5],'color',whichSeparation);
    g(1,1).geom_jitter()
    %g(1,1).stat_bin('geom','line')
    g(1,1).no_legend()
    if thisModel == 2
        g(1,1).set_title('1D within digit direction')
    else
        g(1,1).set_title('2D within digit direction')
    end
    
    g(1,2) = gramm('x',whichSeparation,'y',[PD2y;PD3y;PD4y;PD5y],'color',whichSeparation);
    g(1,2).geom_jitter()
    g(1,2).no_legend()
    
    if thisModel == 2
        g(1,2).set_title('1D between digit direction')
    else
        g(1,2).set_title('2D between digit direction')
    end
    
    
    g.set_names('x',[],'y', 'pRF size')
    g.set_text_options('Font','Helvetica', 'base_size', 12)
    g.set_point_options('base_size',8)
    g.axe_property('XGrid','on','YGrid','on')
    if thisModel == 2
        g.axe_property('YLim',[0 2])
    else
        g.axe_property('YLim',[0 5])
    end
    g.set_order_options('x',0)
    g.draw
    
    g(1,1).update
    g(1,1).stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
    g(1,1).no_legend()
    
    
    g(1,2).update
    g(1,2).stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
    g(1,2).no_legend()
    
    g.draw
    filename = sprintf(['sue_compare_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
else
    clear g
    figure
    
    g = gramm('x',whichSeparation,'y',[PD2;PD3;PD4;PD5],'color',whichSeparation);
    
    g.geom_jitter()
    g.no_legend()
    
    %     if thisModel == 2
    %         g.set_title('1D Within digit model')
    %     elseif thisModel == 3
    %         g.set_title('1D Between digit model')
    %     elseif thisModel == 4
    g.set_title('Unconstrained model')
%    end
    
    g.set_names('x',[],'y', 'pRF size')
    g.set_text_options('Font','Helvetica', 'base_size', 12)
    g.set_point_options('base_size',8)
    g.axe_property('XGrid','on','YGrid','on')
    %g.axe_property('YLim',[0 2])
    g.set_order_options('x',0)
    g.draw
    
    g.update
    g.stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
    g.no_legend()
    
    g.draw
    filename = sprintf(['sue_compare_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
    
    
    
end

%% TRY PDF Version

if strcmpi(sep,'BAA')
    whichSeparation = DIGDEX3; %BA
elseif strcmpi(sep,'BASETIPS')
    whichSeparation = DIGDEX2; %Basetips
end

if thisModel == 1 || thisModel == 2
    clear g
    figure('Position',[100 100 800 400])
    g(1,1) = gramm('x',[PD2;PD3;PD4;PD5],'color',whichSeparation);
    %g(1,1).geom_jitter()
    g(1,1).stat_density()
    %g(1,1).stat_bin('geom','line')
    %g(1,1).no_legend()
    if thisModel == 2
        g(1,1).set_title('1D within digit direction')
    else
        g(1,1).set_title('2D within digit direction')
    end
    
    g(1,2) = gramm('x',[PD2y;PD3y;PD4y;PD5y],'color',whichSeparation);
    g(1,2).stat_density()
    %g(1,2).no_legend()
    
    
    if thisModel == 2
        g(1,2).set_title('1D between digit direction')
    else
        g(1,2).set_title('2D between digit direction')
    end
    
    
    g.set_names('x', 'pRF size') 
    g.set_text_options('Font','Helvetica', 'base_size', 16)
    g.set_point_options('base_size',8)
    g.axe_property('XGrid','on','YGrid','on')
    if thisModel == 1
        g.axe_property('YLim',[0 1.5],'XLim',[0 5])
    else
        g.axe_property('YLim',[0 3.5],'XLim',[0 2.5])
    end
    g.set_order_options('x',0)

    
    g.draw
    %axis square

    filename = sprintf(['sue_compare_' sep '_' datafilenames{thisModel} 'statdensity.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
else
    clear g
    %figure('Position',[100 100 800 400])
    figure
    
    g = gramm('x',[PD2;PD3;PD4;PD5],'color',whichSeparation);
    
    %g.geom_jitter()
    g.stat_density()
    %g.no_legend()
    
    %     if thisModel == 2
    %         g.set_title('1D Within digit model')
    %     elseif thisModel == 3
    %         g.set_title('1D Between digit model')
    %     elseif thisModel == 4
    g.set_title('Unconstrained model')
%    end
    
    g.set_names('x', 'pRF size')
    g.set_text_options('Font','Helvetica', 'base_size', 12)
    g.set_point_options('base_size',8)
    g.axe_property('XGrid','on','YGrid','on')
    %g.axe_property('YLim',[0 2])
    g.set_order_options('x',0)
 
    

    
    g.draw
    filename = sprintf(['sue_compare_' sep '_' datafilenames{thisModel} 'statdensity.pdf'],'%s%s');
    g.export('file_name',filename, ...
        'export_path',...
        '/Users/ppzma/The University of Nottingham/Michael_Sue - pRF/pRF_2021/',...
        'file_type','pdf')
    
    
    
    
end

%% stats



%bastack = [PD2;PD3;PD4;PD5]; % for 2D this is within digit dir
bastack = [PD2y;PD3y;PD4y;PD5y]; % for 2D this is between digit dir
digstack = [D2;D3;D4;D5];
DIGDEXbloop = [repmat({'D2'},nsubs.*4,1); repmat({'D3'},nsubs.*4,1); repmat({'D4'},nsubs.*4,1);repmat({'D5'},nsubs.*4,1)];


if strcmpi(sep,'BAA')
    bawhich = DIGDEX3;
elseif strcmpi(sep,'BASETIPS')
    bawhich = DIGDEX2;
end


%pdstack = [PD2;PD3;PD4;PD5];
%pdwhich = DIGDEX2;

% stats
[p,tbl,stats,terms] = anovan(bastack,{bawhich},'model','linear','varnames',{'BA'}); %interaction?
[resultsBA,~,~,gnamesBA] = multcompare(stats,"Dimension",[1 1]);

tblBA = array2table(resultsBA,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblBA.("Group A")=gnamesBA(tblBA.("Group A"));
tblBA.("Group B")=gnamesBA(tblBA.("Group B"));

writecell(tbl,sprintf(['tbl_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')
writetable(tblBA,sprintf(['mult_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')

%





%% digits

%digstack = [D2;D3;D4;D5];
digstack = [D2y;D3y;D4y;D5y];
[p,tbl,stats,terms] = anovan(digstack,{DIGDEXbloop},'model','linear','varnames',{'Digit'});
[resultsBA,~,~,gnamesBA] = multcompare(stats,"Dimension",[1 1]);

tblBA = array2table(resultsBA,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblBA.("Group A")=gnamesBA(tblBA.("Group A"));
tblBA.("Group B")=gnamesBA(tblBA.("Group B"));

writecell(tbl,sprintf(['digtbl_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')
writetable(tblBA,sprintf(['digmult_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')
















