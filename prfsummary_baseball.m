%% try to refresh prfstat_testingBaseTips % prfstat_testingBAll
% April 2022
% ma
%
%
% this combines prfstat_testingBaseTips and prfstat_testingBAll and
% prfstat_testing_1dwang into one script.
% this is being edited to include 7 BAs, not just 4, but with the option to
% use only 4 PDs.
% this only makes sense if we use a large ROI encompassing precentral
% CG,i.e. underscored _big in datafilenames.
%
% if you want to use the same pool for PD, based on a only post central
% gyrus respones, use prfsummary_baseball_original


% here the beginning is the same as prfsummary2.m
thisorange = [0.85 0.325 0.098];
thisyellow = [0.929 0.694 0.125];
thisgreen = [0.466 0.6740 0.188];
thisblue = [0 0.447 0.741];
thislblue = [0.3010 0.7450 0.9330];
thisred = [0.6350 0.0780 0.1840];
thispurple = [0.4940 0.1840 0.5560];
thiswhite = [1 1 1];
mycmap = [thisorange; thisyellow; thisgreen; thisblue];
mycmap2 = [thisorange; thisyellow; thisgreen; thisblue; thislblue; thisred; thispurple];


mypath='pRF_somato_data/data/';

% range of data
% SUBJECTS
mysubs = {'prf1','prf2','prf3','prf4','prf6','prf7','prf8','prf10'};
% ELIMINATE PRF9 and PRF11 as they don't have consistnet ROIs in
% precentralgyrus
% test case
% mysubs = {'prf1','prf2'};
nsubs = 8;
%howmanyrows = 4;

thisModel = 1;
doBa = 1;
doBasetips = 0;

if doBasetips == 1
    howmanyrows = 4;
elseif doBa == 1
    howmanyrows = 7;
end

% DIGITS / PD locations
digits = 2:5;
baa = 2:howmanyrows+1;
mynbins = 4;
mycutoff = 48;
upperthresh = 15;

datafilenames = {'prf_2ddouble_dd_nov_big', ...
    'prf_1d_nodiag2_nov_big', ...
    'prf_1dt_nodiag2_nov_big', ...
    'prf_sixteen_nodiag2_com_nov_big'};

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

save(sprintf('extraBA_vols_%s_%s.mat',sep,currentModel),'BA_vols_models')


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
                    currentROI = data{iSub}.model(thisModel).roiIndexPD(iLoc, iDigit).idx;
                end
                
                theR2 = currentO(:,myadjr);
                
                [commonVoxels, iInData, iInROI] = intersect(currentIDX, currentROI );
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
    thisModel = 2;
end

%% now plot
% from here, it goes
% Column is a subject
% Rows: D2 D2 D2 D2 D3 etc.
%     :pd2 pd3 pd4 pd4 pd2 etc...

% unpack and average across subjects
for iSub = 1:nsubs
    rfx_sizes = struct2cell(results{iSub}.rfx);
    rfx_sizes = rfx_sizes(:);
    rfx_means(:,iSub) = cellfun(@nanmean,rfx_sizes);
    
    if thisModel == 1
        rfy_sizes = struct2cell(results{iSub}.rfy);
        rfy_sizes = rfy_sizes(:);
        rfy_means(:,iSub) = cellfun(@nanmean,rfy_sizes);
    elseif thisModel == 2 || thisModel == 3
        rfy_sizes = struct2cell(resultsT{iSub}.rfx);
        rfy_sizes = rfy_sizes(:);
        rfy_means(:,iSub) = cellfun(@nanmean,rfy_sizes);
    end
end

if doBasetips == 1
    
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
    end
    
elseif doBa == 1
    
    d2pd2 = rfx_means(1,:);
    d2pd3 = rfx_means(2,:);
    d2pd4 = rfx_means(3,:);
    d2pd5 = rfx_means(4,:);
    d2pd6 = rfx_means(5,:);
    d2pd7 = rfx_means(6,:);
    d2pd8 = rfx_means(7,:);
    
    d3pd2 = rfx_means(8,:);
    d3pd3 = rfx_means(9,:);
    d3pd4 = rfx_means(10,:);
    d3pd5 = rfx_means(11,:);
    d3pd6 = rfx_means(12,:);
    d3pd7 = rfx_means(13,:);
    d3pd8 = rfx_means(14,:);
    
    d4pd2 = rfx_means(15,:);
    d4pd3 = rfx_means(16,:);
    d4pd4 = rfx_means(17,:);
    d4pd5 = rfx_means(18,:);
    d4pd6 = rfx_means(19,:);
    d4pd7 = rfx_means(20,:);
    d4pd8 = rfx_means(21,:);
    
    d5pd2 = rfx_means(22,:);
    d5pd3 = rfx_means(23,:);
    d5pd4 = rfx_means(24,:);
    d5pd5 = rfx_means(25,:);
    d5pd6 = rfx_means(26,:);
    d5pd7 = rfx_means(27,:);
    d5pd8 = rfx_means(28,:);
    
    PD2m = [mean(d2pd2(:)); mean(d3pd2(:)); mean(d4pd2(:)); mean(d5pd2(:))];
    PD3m = [mean(d2pd3(:)); mean(d3pd3(:)); mean(d4pd3(:)); mean(d5pd3(:))];
    PD4m = [mean(d2pd4(:)); mean(d3pd4(:)); mean(d4pd4(:)); mean(d5pd4(:))];
    PD5m = [mean(d2pd5(:)); mean(d3pd5(:)); mean(d4pd5(:)); mean(d5pd5(:))];
    PD6m = [mean(d2pd6(:)); mean(d3pd6(:)); mean(d4pd6(:)); mean(d5pd6(:))];
    PD7m = [mean(d2pd7(:)); mean(d3pd7(:)); mean(d4pd7(:)); mean(d5pd7(:))];
    PD8m = [mean(d2pd8(:)); mean(d3pd8(:)); mean(d4pd8(:)); mean(d5pd8(:))];
    
    PD2mstd = [std(d2pd2(:))./sqrt(length(d2pd2)); std(d3pd2(:))./sqrt(length(d3pd2)); std(d4pd2(:))./sqrt(length(d4pd2)); std(d5pd2(:))./sqrt(length(d5pd2))];
    PD3mstd = [std(d2pd3(:))./sqrt(length(d2pd3)); std(d3pd3(:))./sqrt(length(d3pd3)); std(d4pd3(:))./sqrt(length(d4pd3)); std(d5pd3(:))./sqrt(length(d5pd3))];
    PD4mstd = [std(d2pd4(:))./sqrt(length(d2pd4)); std(d3pd4(:))./sqrt(length(d3pd4)); std(d4pd4(:))./sqrt(length(d4pd4)); std(d5pd4(:))./sqrt(length(d5pd4))];
    PD5mstd = [std(d2pd5(:))./sqrt(length(d2pd5)); std(d3pd5(:))./sqrt(length(d3pd5)); std(d4pd5(:))./sqrt(length(d4pd5)); std(d5pd5(:))./sqrt(length(d5pd5))];
    PD6mstd = [std(d2pd6(:))./sqrt(length(d2pd6)); std(d3pd6(:))./sqrt(length(d3pd6)); std(d4pd6(:))./sqrt(length(d4pd6)); std(d5pd6(:))./sqrt(length(d5pd6))];
    PD7mstd = [std(d2pd7(:))./sqrt(length(d2pd7)); std(d3pd7(:))./sqrt(length(d3pd7)); std(d4pd7(:))./sqrt(length(d4pd7)); std(d5pd7(:))./sqrt(length(d5pd7))];
    PD8mstd = [std(d2pd8(:))./sqrt(length(d2pd8)); std(d3pd8(:))./sqrt(length(d3pd8)); std(d4pd8(:))./sqrt(length(d4pd8)); std(d5pd8(:))./sqrt(length(d5pd8))];
    
    
    PD2 = [d2pd2(:); d3pd2(:); d4pd2(:); d5pd2(:)];
    PD3 = [d2pd3(:); d3pd3(:); d4pd3(:); d5pd3(:)];
    PD4 = [d2pd4(:); d3pd4(:); d4pd4(:); d5pd4(:)];
    PD5 = [d2pd5(:); d3pd5(:); d4pd5(:); d5pd5(:)];
    PD6 = [d2pd6(:); d3pd6(:); d4pd6(:); d5pd6(:)];
    PD7 = [d2pd7(:); d3pd7(:); d4pd7(:); d5pd7(:)];
    PD8 = [d2pd8(:); d3pd8(:); d4pd8(:); d5pd8(:)];
    
    D2 = [d2pd6(:); d2pd7(:); d2pd8(:)];
    D3 = [d3pd6(:); d3pd7(:); d3pd8(:)];
    D4 = [d4pd6(:); d4pd7(:); d4pd8(:)];
    D5 = [d5pd6(:); d5pd7(:); d5pd8(:)];
    
    
    if thisModel == 1 || thisModel == 2 || thisModel == 3
        d2pd2y = rfy_means(1,:);
        d2pd3y = rfy_means(2,:);
        d2pd4y = rfy_means(3,:);
        d2pd5y = rfy_means(4,:);
        d2pd6y = rfy_means(5,:);
        d2pd7y = rfy_means(6,:);
        d2pd8y = rfy_means(7,:);
        
        d3pd2y = rfy_means(8,:);
        d3pd3y = rfy_means(9,:);
        d3pd4y = rfy_means(10,:);
        d3pd5y = rfy_means(11,:);
        d3pd6y = rfy_means(12,:);
        d3pd7y = rfy_means(13,:);
        d3pd8y = rfy_means(14,:);
        
        d4pd2y = rfy_means(15,:);
        d4pd3y = rfy_means(16,:);
        d4pd4y = rfy_means(17,:);
        d4pd5y = rfy_means(18,:);
        d4pd6y = rfy_means(19,:);
        d4pd7y = rfy_means(20,:);
        d4pd8y = rfy_means(21,:);
        
        d5pd2y = rfy_means(22,:);
        d5pd3y = rfy_means(23,:);
        d5pd4y = rfy_means(24,:);
        d5pd5y = rfy_means(25,:);
        d5pd6y = rfy_means(26,:);
        d5pd7y = rfy_means(27,:);
        d5pd8y = rfy_means(28,:);
        
        PD2ym = [mean(d2pd2y(:)); mean(d3pd2y(:)); mean(d4pd2y(:)); mean(d5pd2y(:))];
        PD3ym = [mean(d2pd3y(:)); mean(d3pd3y(:)); mean(d4pd3y(:)); mean(d5pd3y(:))];
        PD4ym = [mean(d2pd4y(:)); mean(d3pd4y(:)); mean(d4pd4y(:)); mean(d5pd4y(:))];
        PD5ym = [mean(d2pd5y(:)); mean(d3pd5y(:)); mean(d4pd5y(:)); mean(d5pd5y(:))];
        PD6ym = [mean(d2pd6y(:)); mean(d3pd6y(:)); mean(d4pd6y(:)); mean(d5pd6y(:))];
        PD7ym = [mean(d2pd7y(:)); mean(d3pd7y(:)); mean(d4pd7y(:)); mean(d5pd7y(:))];
        PD8ym = [mean(d2pd8y(:)); mean(d3pd8y(:)); mean(d4pd8y(:)); mean(d5pd8y(:))];
        
        PD2ymstd = [std(d2pd2y(:))./sqrt(length(d2pd2y)); std(d3pd2y(:))./sqrt(length(d3pd2y)); std(d4pd2y(:))./sqrt(length(d4pd2y)); std(d5pd2y(:))./sqrt(length(d5pd2y))];
        PD3ymstd = [std(d2pd3y(:))./sqrt(length(d2pd3y)); std(d3pd3y(:))./sqrt(length(d3pd3y)); std(d4pd3y(:))./sqrt(length(d4pd3y)); std(d5pd3y(:))./sqrt(length(d5pd3y))];
        PD4ymstd = [std(d2pd4y(:))./sqrt(length(d2pd4y)); std(d3pd4y(:))./sqrt(length(d3pd4y)); std(d4pd4y(:))./sqrt(length(d4pd4y)); std(d5pd4y(:))./sqrt(length(d5pd4y))];
        PD5ymstd = [std(d2pd5y(:))./sqrt(length(d2pd5y)); std(d3pd5y(:))./sqrt(length(d3pd5y)); std(d4pd5y(:))./sqrt(length(d4pd5y)); std(d5pd5y(:))./sqrt(length(d5pd5y))];
        PD6ymstd = [std(d2pd6y(:))./sqrt(length(d2pd6y)); std(d3pd6y(:))./sqrt(length(d3pd6y)); std(d4pd6y(:))./sqrt(length(d4pd6y)); std(d5pd6y(:))./sqrt(length(d5pd6y))];
        PD7ymstd = [std(d2pd7y(:))./sqrt(length(d2pd7y)); std(d3pd7y(:))./sqrt(length(d3pd7y)); std(d4pd7y(:))./sqrt(length(d4pd7y)); std(d5pd7y(:))./sqrt(length(d5pd7y))];
        PD8ymstd = [std(d2pd8y(:))./sqrt(length(d2pd8y)); std(d3pd8y(:))./sqrt(length(d3pd8y)); std(d4pd8y(:))./sqrt(length(d4pd8y)); std(d5pd8y(:))./sqrt(length(d5pd8y))];
        
        
        PD2y = [d2pd2y(:); d3pd2y(:); d4pd2y(:); d5pd2y(:)];
        PD3y = [d2pd3y(:); d3pd3y(:); d4pd3y(:); d5pd3y(:)];
        PD4y = [d2pd4y(:); d3pd4y(:); d4pd4y(:); d5pd4y(:)];
        PD5y = [d2pd5y(:); d3pd5y(:); d4pd5y(:); d5pd5y(:)];
        PD6y = [d2pd6y(:); d3pd6y(:); d4pd6y(:); d5pd6y(:)];
        PD7y = [d2pd7y(:); d3pd7y(:); d4pd7y(:); d5pd7y(:)];
        PD8y = [d2pd8y(:); d3pd8y(:); d4pd8y(:); d5pd8y(:)];
        
        
        D2y = [d2pd6y(:); d2pd7y(:); d2pd8y(:)];
        D3y = [d3pd6y(:); d3pd7y(:); d3pd8y(:)];
        D4y = [d4pd6y(:); d4pd7y(:); d4pd8y(:)];
        D5y = [d5pd6y(:); d5pd7y(:); d5pd8y(:)];
    
    end
    
end

%
DIGDEX = [repmat({'D2'},nsubs,1); repmat({'D3'},nsubs,1); repmat({'D4'},nsubs,1);repmat({'D5'},nsubs,1)];
DIGDEX2 = [repmat({'PD2'},nsubs.*4,1); repmat({'PD3'},nsubs.*4,1); repmat({'PD4'},nsubs.*4,1);repmat({'PD5'},nsubs.*4,1)];
DIGDEX3 = [repmat({'BA3a'},nsubs.*4,1); repmat({'BA3b'},nsubs.*4,1); repmat({'BA1'},nsubs.*4,1);repmat({'BA2'},nsubs.*4,1)];
DIGDEX33 = [repmat({'BA3a'},nsubs.*4,1); repmat({'BA3b'},nsubs.*4,1); repmat({'BA1'},nsubs.*4,1);repmat({'BA2'},nsubs.*4,1);...
    repmat({'BA4a'},nsubs.*4,1);repmat({'BA4p'},nsubs.*4,1);repmat({'BA6'},nsubs.*4,1)];

%% explode plot
if thisModel == 1 || thisModel == 2 || thisModel == 3
    clear g
    figure('Position',[100 100 1400 400])
    g(1,1) = gramm('x',PD2y,'y',PD2,'color',DIGDEX);
    g(1,1).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,1).set_title('PD2')
    else
        g(1,1).set_title('BA3a')
    end
    g(1,2) = gramm('x',PD3y,'y',PD3,'color',DIGDEX);
    g(1,2).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,2).set_title('PD3')
    else
        g(1,2).set_title('BA3b')
    end
    g(1,3) = gramm('x',PD4y,'y',PD4,'color',DIGDEX);
    g(1,3).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,3).set_title('PD4')
    else
        g(1,3).set_title('BA1')
    end
    g(1,4) = gramm('x',PD5y,'y',PD5,'color',DIGDEX);
    g(1,4).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,4).set_title('PD5')
    else
        g(1,4).set_title('BA2')
    end
    
    if doBa
        g(1,5) = gramm('x',PD6y,'y',PD6,'color',DIGDEX);
        g(1,5).geom_jitter()
        g(1,5).set_title('BA4a')
        g(1,6) = gramm('x',PD7y,'y',PD7,'color',DIGDEX);
        g(1,6).geom_jitter()
        g(1,6).set_title('BA4p')
        g(1,7) = gramm('x',PD8y,'y',PD8,'color',DIGDEX);
        g(1,7).geom_jitter()
        g(1,7).set_title('BA6')
        
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
    
    if doBa
        g(1,5).update('x',x,'y',y)
        g(1,5).geom_line
        g(1,6).update('x',x,'y',y)
        g(1,6).geom_line
        g(1,7).update('x',x,'y',y)
        g(1,7).geom_line
    end
    
    g.set_color_options('map',[0.5 0.5 0.5])
    g.draw()
    filename = sprintf(['explode_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
%     g.export('file_name',filename, ...
%         'export_path',...
%         '',...
%         'file_type','pdf')
    
    
    
else
    %DIGDEX3 = [{'PD2'},{'PD3'},{'PD4'},{'PD5'}];
    %DIGDEX4 = [repmat(DIGDEX3,4,1)];
    
    clear g
    figure('Position',[100 100 1400 400])
    g(1,1) = gramm('x',DIGDEX,'y',PD2,'color',DIGDEX);
    g(1,1).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,1).set_title('PD2')
    else
        g(1,1).set_title('BA3a')
    end
    g(1,2) = gramm('x',DIGDEX,'y',PD3,'color',DIGDEX);
    g(1,2).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,2).set_title('PD3')
    else
        g(1,2).set_title('BA3b')
    end
    g(1,3) = gramm('x',DIGDEX,'y',PD4,'color',DIGDEX);
    g(1,3).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,3).set_title('PD4')
    else
        g(1,3).set_title('BA1')
    end
    g(1,4) = gramm('x',DIGDEX,'y',PD5,'color',DIGDEX);
    g(1,4).geom_jitter()
    if strcmpi(sep,'BASETIPS')
        g(1,4).set_title('PD5')
    else
        g(1,4).set_title('BA2')
    end
    
    if doBa
        g(1,5) = gramm('x',DIGDEX,'y',PD6,'color',DIGDEX);
        g(1,5).geom_jitter()
        g(1,5).set_title('BA4a')
        g(1,6) = gramm('x',DIGDEX,'y',PD7,'color',DIGDEX);
        g(1,6).geom_jitter()
        g(1,6).set_title('BA4p')
        g(1,7) = gramm('x',DIGDEX,'y',PD8,'color',DIGDEX);
        g(1,7).geom_jitter()
        g(1,7).set_title('BA6')
        
    end
    
    
    g.set_names('x','Digit','y','Moment of Inertia')
    g.set_text_options('Font','Helvetica', 'base_size', 16)
    g.set_point_options('base_size',10)
    g.axe_property('XGrid','on','YGrid','on','YLim',[5 8],'DataAspectRatio',[1 1 1])
    g.set_order_options('x',0)
    g.set_color_options('map',mycmap)
    g.draw
    filename = sprintf(['explode_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
%     g.export('file_name',filename, ...
%         'export_path',...
%         '',...
%         'file_type','pdf')
    
    
    
end

%% now collapse above plot
aa = 2.4; %axis limits 1D is 2.4, 2D Gaussian is 3.5
% stat_summary('geom',{'point' 'errorbar'},'dodge',0.3,'width',0.5);
DIGDEXm = [repmat({'D2'},1,1); repmat({'D3'},1,1); repmat({'D4'},1,1);repmat({'D5'},1,1)];
if strcmpi(sep,'BASETIPS')
    DIGDEX2m = [repmat({'PD2'},4,1); repmat({'PD3'},4,1); repmat({'PD4'},4,1);repmat({'PD5'},4,1)];
else
    DIGDEX2m = [repmat({'BA3a'},4,1); repmat({'BA3b'},4,1); repmat({'BA1'},4,1);repmat({'BA2'},4,1)];
    DIGDEX3m = [repmat({'BA3a'},4,1); repmat({'BA3b'},4,1); repmat({'BA1'},4,1);repmat({'BA2'},4,1);...
        repmat({'BA4a'},4,1);repmat({'BA4p'},4,1);repmat({'BA6'},4,1)];
    
end

mycolor1 = [thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue];

% mycolor2 = [thisorange;thisyellow;thisgreen;thisblue;thislblue;thisred;thispurple;...
%     thisorange;thisyellow;thisgreen;thisblue;thislblue;thisred;thispurple;...
%     thisorange;thisyellow;thisgreen;thisblue;thislblue;thisred;thispurple;...
%     thisorange;thisyellow;thisgreen;thisblue;thislblue;thisred;thispurple;];

mycolor2 = [thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
    thisorange;thisyellow;thisgreen;thisblue;...
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

mymarker2 = ['s';'s';'s';'s';...
    '^';'^';'^';'^';...
    'o';'o';'o';'o';...
    'p';'p';'p';'p';...
    'v';'v';'v';'v';...
    '>';'>';'>';'>';...
    '<';'<';'<';'<'];


if thisModel == 1 || thisModel == 2 || thisModel == 3
    clear g
    %figure('Position',[100 100 1400 400])
    figure
    if doBasetips == 1
        g = gramm('x',[PD2ym;PD3ym;PD4ym;PD5ym],'y',[PD2m;PD3m;PD4m;PD5m],'color',[DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm],'marker',DIGDEX2m);
    elseif doBa == 1
        g = gramm('x',[PD2ym;PD3ym;PD4ym;PD5ym;PD6ym;PD7ym;PD8ym],'y',[PD2m;PD3m;PD4m;PD5m;PD6m;PD7m;PD8m],'color',[DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm],'marker',DIGDEX3m);
    end
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
%     g.export('file_name',filename, ...
%         'export_path',...
%         '',...
%         'file_type','pdf')
    
    % wang plot
    if doBa == 1
        DYstack = [PD2ym;PD3ym;PD4ym;PD5ym;PD6ym;PD7ym;PD8ym];
        DXstack = [PD2m;PD3m;PD4m;PD5m;PD6m;PD7m;PD8m];
        DYstack_std = [PD2ymstd;PD3ymstd;PD4ymstd;PD5ymstd;PD6ymstd;PD7ymstd;PD8ymstd];
        DXstack_std = [PD2mstd;PD3mstd;PD4mstd;PD5mstd;PD6mstd;PD7mstd;PD8mstd];
    elseif doBasetips == 1
        DYstack = [PD2ym;PD3ym;PD4ym;PD5ym];
        DXstack = [PD2m;PD3m;PD4m;PD5m];
        DYstack_std = [PD2ymstd;PD3ymstd;PD4ymstd;PD5ymstd];
        DXstack_std = [PD2mstd;PD3mstd;PD4mstd;PD5mstd];
    end
    
    figure
    if doBa == 1
        for kk = 17:28 % for just 4a,4p and 6, do 17:28
            e = errorbar(DYstack(kk),DXstack(kk), DYstack_std(kk), DXstack_std(kk),...
                'both','o','MarkerEdgeColor', mycolor2(kk,:),'MarkerSize',10,...
                'MarkerFaceColor',mycolor2(kk,:),'Marker',mymarker2(kk));
            e.Color = mycolor2(kk,:);
            %e.Color = [0.6 0.6 0.6];
            hold on
        end
    elseif doBasetips == 1
        for kk = 1:16
            e = errorbar(DYstack(kk),DXstack(kk), DYstack_std(kk), DXstack_std(kk),...
                'both','o','MarkerEdgeColor', mycolor1(kk,:),'MarkerSize',10,...
                'MarkerFaceColor',mycolor1(kk,:),'Marker',mymarker(kk));
            e.Color = mycolor1(kk,:);
            hold on
        end
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
    
    filename = sprintf(['_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    %print('-dpdf', filename)
    
    [H,P,CI,STATS] = ttest2(DYstack, DXstack);

    
elseif thisModel == 4
    clear g
    %figure('Position',[100 100 1400 400])
    figure
    if doBasetips == 1
        g = gramm('x',DIGDEX2m,'y',[PD2m;PD3m;PD4m;PD5m],'color',[DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm]);
    elseif doBa == 1
        g = gramm('x',DIGDEX3m,'y',[PD2m;PD3m;PD4m;PD5m;PD6m;PD7m;PD8m],'color',[DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm;DIGDEXm]);
    end
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
%     g.export('file_name',filename, ...
%         'export_path',...
%         '',...
%         'file_type','pdf')
    
    
    if doBa == 1
        DXstack = [PD2m;PD3m;PD4m;PD5m;PD6m;PD7m;PD8m];
        DXstack_std = [PD2mstd;PD3mstd;PD4mstd;PD5mstd;PD6mstd;PD7mstd;PD8mstd];
        bigcolor2 = [1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5,3.2,3.3,3.4,3.5,4.2,4.3,4.4,4.5,5.2,5.3,5.4,5.5,6.2,6.3,6.4,6.5,7.2,7.3,7.4,7.5];
        bigcolor3 = [1.35 2.35 3.35 4.35 5.35 6.35 7.35];
    elseif doBasetips == 1
        DXstack = [PD2m;PD3m;PD4m;PD5m];
        DXstack_std = [PD2mstd;PD3mstd;PD4mstd;PD5mstd];
        bigcolor2 = [1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5,3.2,3.3,3.4,3.5,4.2,4.3,4.4,4.5];
        bigcolor3 = [1.35 2.35 3.35 4.35];
    end
    
    figure()
    hold on
    if doBa == 1
        for kk = 1:28
            e = errorbar(bigcolor2(kk),DXstack(kk), DXstack_std(kk),...
                'o','MarkerEdgeColor', mycolor2(kk,:),'MarkerSize',10,...
                'MarkerFaceColor',mycolor2(kk,:));
            e.Color = mycolor2(kk,:);
            hold on
        end
        xticks([1.35 2.35 3.35 4.35 5.35 6.35 7.35])
        xticklabels({'BA3a','BA3b','BA1','BA2','BA4a','BA4p','BA6'})
        xlim([1 7.5])
    elseif doBasetips == 1
        for kk = 1:16
            e = errorbar(bigcolor2(kk),DXstack(kk), DXstack_std(kk),...
                'o','MarkerEdgeColor', mycolor1(kk,:),'MarkerSize',10,...
                'MarkerFaceColor',mycolor1(kk,:));
            e.Color = mycolor1(kk,:);
            hold on
        end
        xticks([1.35 2.35 3.35 4.35])
        xticklabels({'PD2','PD3','PD4','PD5'})
        xlim([1 4.5])
    end
    
    
    
    
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
    filename = sprintf(['_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
    %print('-dpdf', filename)
    
end


%% try all subjects wang plot
%
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

%% now  plot

if strcmpi(sep,'BAA')
    whichSeparation = DIGDEX33; %BA
    howlong = 1400;
elseif strcmpi(sep,'BASETIPS')
    whichSeparation = DIGDEX2; %Basetips
    howlong = 500;
end



if thisModel == 1 || thisModel == 2
    clear g
    figure('Position',[100 100 howlong 400])
    if doBasetips == 1
        g(1,1) = gramm('x',DIGDEX2,'y',[PD2;PD3;PD4;PD5],'color',DIGDEX2);
        g(1,1).geom_jitter()
        g(1,1).no_legend()
        
        g(1,2) = gramm('x',DIGDEX2,'y',[PD2y;PD3y;PD4y;PD5y],'color',DIGDEX2);
        g(1,2).geom_jitter()
        g(1,2).no_legend()
        
        if thisModel == 2
            g(1,1).set_title('1D within digit direction')
            g(1,2).set_title('1D between digit direction')  
        elseif thisModel == 1
            g(1,1).set_title('2D within digit direction')
            g(1,2).set_title('2D between digit direction')
        end
    
            
    else
        g(1,1) = gramm('x',DIGDEX33,'y',[PD2;PD3;PD4;PD5;PD6;PD7;PD8],'color',DIGDEX33);
        g(1,1).geom_jitter()
        g(1,1).no_legend()
        
        g(1,2) = gramm('x',DIGDEX33,'y',[PD2y;PD3y;PD4y;PD5y;PD6y;PD7y;PD8y],'color',DIGDEX33);
        g(1,2).geom_jitter()
        g(1,2).no_legend()
        
        if thisModel == 2
            g(1,1).set_title('1D within digit direction')
            g(1,2).set_title('1D between digit direction')
        elseif thisModel == 1
            g(1,1).set_title('2D within digit direction')
            g(1,2).set_title('2D between digit direction')
        end
        
                
    end
    
    
    g.set_names('x',[],'y', 'pRF size')
    g.set_text_options('Font','Helvetica', 'base_size', 12)
    g.set_point_options('base_size',8)
    g.axe_property('XGrid','on','YGrid','on')
    %g.axe_property('YLim',[0 6])
    if thisModel == 2
        g.axe_property('YLim',[0 2])
    else
        g.axe_property('YLim',[0 6])
    end

    g.set_order_options('x',0)
    g.draw
    
    g(1,1).update
    g(1,1).stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
    g(1,1).no_legend()
    %g(1,1).geom_line()
    
    g(1,2).update
    g(1,2).stat_boxplot('width', 2, 'dodge', 0.5, 'alpha', 0.5, 'linewidth', 2, 'drawoutlier',0)
    g(1,2).no_legend()
    
    g.draw
    filename = sprintf(['sue_compare_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
%     g.export('file_name',filename, ...
%         'export_path',...
%         '',...
%         'file_type','pdf')
    
else
    clear g
    figure('Position',[100 100 howlong 400])
    if doBa == 1
        if thisModel == 2 || thisModel == 3
            g = gramm('x',DIGDEX33,'y',[PD2y;PD3y;PD4y;PD5y;PD6y;PD7y;PD8y],'color',DIGDEX33);
        else
            g = gramm('x',DIGDEX33,'y',[PD2;PD3;PD4;PD5;PD6;PD7;PD8],'color',DIGDEX33);
        end
    elseif doBasetips == 1
        if thisModel == 2 || thisModel == 3
            g = gramm('x',DIGDEX2,'y',[PD2y;PD3y;PD4y;PD5y],'color',DIGDEX2);
        else
            g = gramm('x',DIGDEX2,'y',[PD2;PD3;PD4;PD5],'color',DIGDEX2);
        end
    end
    g.geom_jitter()
    g.no_legend()
    
    if thisModel == 2
        g.set_title('1D Within digit model')
    elseif thisModel == 3
        g.set_title('1D Between digit model')
    elseif thisModel == 4
        g.set_title('Unconstrained model')
    end
    
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
    filename = sprintf(['compare_' sep '_' datafilenames{thisModel} '.pdf'],'%s%s');
%     g.export('file_name',filename, ...
%         'export_path',...
%         '',...
%         'file_type','pdf')
%     
    
    
    
end




%%
%% stats



bastack = [PD6;PD7;PD8];

% if strcmpi(sep,'BAA')
%     bawhich = DIGDEX3;
% elseif strcmpi(sep,'BASETIPS')
%     bawhich = DIGDEX2;
% end
digstack = [D2;D3;D4;D5];

DIGDEXbloop = [repmat({'D2'},length(D2),1); repmat({'D3'},length(D3),1); repmat({'D4'},length(D4),1); repmat({'D5'},length(D5),1)];

%DIGDEXbloop = [repmat({'D2'},nsubs.*4,1); repmat({'D3'},nsubs.*4,1); repmat({'D4'},nsubs.*4,1);repmat({'D5'},nsubs.*4,1)];


bawhich = DIGDEX33(129:end);

%pdstack = [PD2;PD3;PD4;PD5];
%pdwhich = DIGDEX2;

% stats
[p,tbl,stats,terms] = anovan(bastack,{bawhich},'model','interaction','varnames',{'BA'});
[resultsBA,~,~,gnamesBA] = multcompare(stats,"Dimension",[1 1]);

tblBA = array2table(resultsBA,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblBA.("Group A")=gnamesBA(tblBA.("Group A"));
tblBA.("Group B")=gnamesBA(tblBA.("Group B"));

%writecell(tbl,sprintf(['tbl_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')
%writetable(tblBA,sprintf(['mult_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')




% digits
[p,tbl,stats,terms] = anovan(digstack,{DIGDEXbloop},'model','interaction','varnames',{'Digit'});
[resultsBA,~,~,gnamesBA] = multcompare(stats,"Dimension",[1 1]);

tblBA = array2table(resultsBA,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblBA.("Group A")=gnamesBA(tblBA.("Group A"));
tblBA.("Group B")=gnamesBA(tblBA.("Group B"));

%writecell(tbl,sprintf(['digtbl_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')
%writetable(tblBA,sprintf(['digmult_' sep '_' datafilenames{thisModel},'%s%s']),'FileType','spreadsheet')






