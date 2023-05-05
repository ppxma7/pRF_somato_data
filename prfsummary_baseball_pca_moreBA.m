
mypath='/Volumes/nemosine/prfsomato_fits';
%mypath = '/Users/lpzds1/Library/CloudStorage/OneDrive-SharedLibraries-TheUniversityofNottingham/Michael_Sue - pRF_paper_draft/prfsomato_fits';

% range of data
% SUBJECTS
mysubs = {'prf1','prf2','prf3','prf4','prf6','prf7','prf8','prf10'};
% test case
% mysubs = {'prf1','prf2'};
nsubs = 8;

% DIGITS / PD locations
digits = 2:5;
locs = 2:8;
mynbins = 4;
%mycutoff = 48;
upperthresh = 15;

datafilenames = {'prf_2ddouble_dd_nov_big', ...
    'prf_1d_nodiag2_nov_big', ...
    'prf_1dt_nodiag2_nov_big', ...
    'prf_sixteen_nodiag2_com_nov_big'};

nModels = numel(datafilenames);

data = {}; % empty struct
results = {};

for iSub = 1:nsubs
    
    currentSub = mysubs{iSub};
    
    tic()
    % for each data file
    for iModel = 1:nModels
        currentModel = datafilenames{iModel};
        
        % Read in all BA ROIs and indices
        for iLoc = 1:7
            currentBaa = locs(iLoc);
            
            for iDigit = 1:4
                %currentBaa = baa(iBaa);
                currentDigit = digits(iDigit);
                % just pick one model (ROI indeces are all the same
                % load([mypath mysubs{iSub} '/baa/' 'D' num2str(iDig) 'pd' num2str(iBaa) '_16.mat'],'prf_overlays');
                currentROIFilename = [mypath filesep() currentSub filesep() '/baa/' sprintf('D%dpd%d_16.mat', currentDigit, currentBaa)];
                data{iSub}.roiIndex(iLoc, iDigit) = load(currentROIFilename, 'idx');
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
        
        if iModel == 4  %2DGaussian
            
            % fudge factor as maps are flipped from pRFsomatoFit output for
            % unconstrained for some reason
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
            %keyboard
            % there is an abs situation here in pRFsomatoFit so they have
            % to be flipped
            pd_2d_1dex = thisPD_conv==1;
            pd_2d_2dex = thisPD_conv==2;
            pd_2d_3dex = thisPD_conv==3;
            pd_2d_4dex = thisPD_conv==4;
            thisPD_conv(pd_2d_1dex) = (thisPD_conv(pd_2d_1dex).*-1)+5;
            thisPD_conv(pd_2d_2dex) = (thisPD_conv(pd_2d_2dex).*-1)+5;
            thisPD_conv(pd_2d_3dex) = (thisPD_conv(pd_2d_3dex).*-1)+5;
            thisPD_conv(pd_2d_4dex) = (thisPD_conv(pd_2d_4dex).*-1)+5;
            
            data{iSub}.model(iModel).convPD = thisPD_conv;
            
        elseif iModel == 1
%             pd_2d_1dex = thisPD_conv==1;
%             pd_2d_2dex = thisPD_conv==2;
%             pd_2d_3dex = thisPD_conv==3;
%             pd_2d_4dex = thisPD_conv==4;
%             thisPD_conv(pd_2d_1dex) = 4;
%             thisPD_conv(pd_2d_2dex) = 3;
%             thisPD_conv(pd_2d_3dex) = 2;
%             thisPD_conv(pd_2d_4dex) = 1;
            
            % 2D are not flipped but the overlays are wrong
            
            data{iSub}.model(iModel).convPD = thisDig_conv;
            data{iSub}.model(iModel).convDig = thisPD_conv;

            
        else
            
            data{iSub}.model(iModel).convPD = thisPD_conv;
        end
        
        % now we can overlap and index to get PD ROIs.
        % this is not like BA, which all have same ROIs for all models
        % for each model the PD ROIs will be different, as they are
        % dependent on where the values lie.
        for iPD = 1:7
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
thisModel = 3;
doBa = 1;
doBasetips = 0;
if doBa == 1
    sep = 'BAA';
    howmanyrows = 7;
elseif doBasetips == 1
    sep = 'BASETIPS';
    howmanyrows = 4;
end


aa = 4; %axis limits

doPCA = 1;
doMean = 0;

%iSub = 10;
%iDigit = 1;
%iBaa = 1;
%iModel = 4; % 1 is 2D Gaussian. 2 is 1D. 3 is 1DT. 4 is unconstrained.

stims = 4;
nSteps = 3; %fudge to make sure we have a 4x4 grid from each model

mycutoff = 48; % how many voxels to include in the rank
mydivs = divisor(mycutoff); % for figure
myrangelen = [(length(mydivs)./2) (length(mydivs)./2)+1];
myrange = [mydivs(myrangelen(1)) mydivs(myrangelen(2))];

whichComponent = 1;
spitout = 1; %display images



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

for iSub = 1:nsubs
    myPCA = struct;
   
    figure('Position',[100 100 1600 1600])
    tiledlayout(howmanyrows,4)
    
    for iLoc = 1:howmanyrows
        for iDigit = 1:4
            
            currentSub = mysubs{iSub};
            currentDigit = digits(iDigit);
            currentBaa = locs(iLoc);
            currentModel = datafilenames{thisModel};
            currentIDX = data{iSub}.model(iModel).data.idx_limited;
            % current data
            currentD = data{iSub}.model(thisModel).data.dd;
            
            %currentO = data{iSub}.model(thisModel).data.prf_overlays_limited;
%            [commonVoxels, iInData, iInROI] = intersect(currentIDX, currentROI );
            % ROI - subset that's in the particular ROI -- same as in linearCoords
            % here is where we edit this for either BAA or PD
            linearCoords = currentD.linearCoords;
            if doBa == 1
                %currentROI = data{iSub}.roiIndex(iDigit, iBaa).idx;
                %currentO = data{iSub}.model(iModel).data.prf_overlays;
                
                currentROI = data{iSub}.roiIndex(iLoc, iDigit).idx;
                [commonVoxels, iInData, iInROI] = intersect(linearCoords, currentROI );
                theR2 = currentD.r2;
            elseif doBasetips == 1
                %keyboard
                % need to create the idx here
                currentROI = data{iSub}.model(thisModel).roiIndexPD(iLoc, iDigit).idx;
                %[commonVoxels, iInData, iInROI] = intersect(currentIDX, currentROI );
                [commonVoxels, iInData, iInROI] = intersect(linearCoords, currentROI );
                %theR2 = currentO(:,myadjr);
                theR2 = currentD.r2;
            end
            
            %theR2 = currentO(:,myadjr);
            
            extractedParams = currentD.params;%(:, ~invalidCoords);
            
            % all DATA / params
            %extractedParams = currentD.params(:, ~invalidCoords);
            
            %[commonVoxels, iInData, iInROI] = intersect(currentIDX, currentROI );
            roiR2 = theR2(iInData);
            roiR2 = roiR2(:);
            
            [topN, topNidx] = maxk(roiR2, mycutoff); % to be consistent with PCA later
            
            matchingLinearIdx = iInData(topNidx);
            %matchingLinearIdx = iInData;
            % okay we don't care about pRFs here, only just getting the
            % indices for each ROI, then we can index into pRF size values
            % and plot these.
            
            % can we save out some TCs for later plotting?
            bestTCs = currentD.mytSeries(:,matchingLinearIdx);
            bestModResp = currentD.mymodelResp(:,matchingLinearIdx);
            bestr2s = currentD.r2(matchingLinearIdx);
            
            theBestR2 = bestr2s;
            theBestTC = bestTCs;%bestTCs(:,1); % get top ranked one
            theBestModResp = bestModResp; %bestModResp(:,1);
            myBestTC.tc{iSub}{iLoc,iDigit} = theBestTC;
            myBestTC.modresp{iSub}{iLoc,iDigit} = theBestModResp;
            myBestTC.bestr2{iSub}{iLoc,iDigit} = theBestR2;
            myBestTC.MLI{iSub}{iLoc,iDigit} = matchingLinearIdx;
            
            if thisModel == 1 %2D Gaussian
                [thisrf,B] = model_2d_gauss(extractedParams(:, matchingLinearIdx) ,stims, nSteps, 1); % extra flag for 1 = smooth

            elseif thisModel == 2 %1D
                [thisrf,B] = model_1d_gauss(extractedParams(:, matchingLinearIdx) ,stims, nSteps+1, 0,1);
            elseif thisModel == 3 %1D T
                [thisrf,B] = model_1d_gauss(extractedParams(:, matchingLinearIdx) ,stims, nSteps+1, 1,1); % extra flag for transpose
            elseif thisModel == 4 % Uncon
                thisrf = model_uncon(extractedParams(:, matchingLinearIdx) ,stims, nSteps, 1);% extra flag for 1 = smooth
            end

            ogcopyrf = thisrf;
            thisrf = B;

            if doPCA == 1
                
                X = reshape(thisrf, [size(thisrf,1).*size(thisrf,1) size(thisrf,3)]);
                % need to flip it as ROWS = Observations, and COLUMNS = Variables;
                % so e.g. we have 48 observations (voxels), and 16 variables (each pRF
                % shape), for one ROI, e.g. D2, PD2
                Xt = transpose(X);
                [coeff, score, latent, tsquared, explained] = pca(Xt);
                
                if isempty(coeff) % error checking in case there's only 1 matching IDX
                    coeff = Xt;
                    coeff = coeff(:);
                    explained = 0;
                end
                %check
                if size(coeff,2)==3 || size(coeff,2)>3
                    % each column of coeff is a principal component in descending order of
                    % component variance
                    topPC = coeff(:,1); % first principal component
                    topPC_rf = reshape(topPC,size(thisrf,1),size(thisrf,2));
                    
                    topPC2 = coeff(:,2); % second PC
                    topPC3 = coeff(:,3); % third PC
                    
                    topPC_rf2 = reshape(topPC2,size(thisrf,1),size(thisrf,2));
                    topPC_rf3 = reshape(topPC3,size(thisrf,1),size(thisrf,2));
                    
                else
                    topPC = coeff(:,1); % first principal component
                    topPC_rf = reshape(topPC,size(thisrf,1),size(thisrf,2));
                    
                end
                % need to label these
                % making a cell array of structs
                myPCA.group{iLoc,iDigit} = coeff;
                
                myPCA.explained{iLoc,iDigit} = explained;
                
                nexttile
                if whichComponent == 1
                    imagesc(topPC_rf);
                elseif whichComponent == 2
                    imagesc(topPC_rf2);
                elseif whichComponent == 3
                    imagesc(topPC_rf3);
                end
                
         
                
                
                axis square
                colormap parula
                t_ = text(-2,-2,sprintf('sub:%s, dig:%d, pd%d, m:%s pca:%d', currentSub, currentDigit, ...
                    currentBaa, currentModel, whichComponent));
                set(t_, 'color', [0 0 0 ], ...
                    'fontsize', 12, 'fontname', 'helvetica', 'interpreter', 'none', ...
                    'verticalalignment', 'top');
                
                text(-2,-3,[num2str(explained(whichComponent)) ' % explained'],'Color','black','FontSize',12)
                
                
                myRF.group{iSub}{iLoc,iDigit} = thisrf;
                
                %means
            elseif doMean == 1
                X = mean(thisrf,3);
                myMEAN.group{iLoc,iDigit} = X;
                nexttile %(iDigit)
                imagesc(X);
                axis square
                colormap parula
                t_ = text(0,0,sprintf('sub: %s, digit: %d, pd %d, m: %s', currentSub, currentDigit, ...
                    currentBaa, currentModel));
                set(t_, 'color', [0 0 0 ], ...
                    'fontsize', 12, 'fontname', 'helvetica', 'interpreter', 'none', ...
                    'verticalalignment', 'top');
                text(1,1,['med r2 ' num2str(median(roiR2))],'Color','black','FontSize',14)
            end
            
                
 
        end
    end
    
    if doPCA == 1
        myGRANDpca.group{iSub} = myPCA.group;
        myGRANDpca.explained{iSub} = myPCA.explained;
        
        print(gcf, sprintf('PCA_baa%d%s_%spca%d%s.png', currentBaa,...
            currentModel, currentSub, whichComponent, sep), '-r300', '-dpng')
        
    elseif doMean == 1
        myGRANDmean.group{iSub} = myMEAN.group;
        print(gcf, sprintf('MEAN_baa%d%s_%spca%d%s.png', currentBaa,...
            currentModel, currentSub, whichComponent, sep), '-r300', '-dpng')
    end
    
end



if doMean == 1
    savefile = ['/Volumes/nemosine/prfsomato_fits/' currentModel sep '_myMEAN_big.mat'];
    save(savefile,'myGRANDmean');
    
elseif doPCA == 1
    savefile = ['/Volumes/nemosine/prfsomato_fits/' currentModel sep '_myPCA_big.mat'];
    save(savefile,'myGRANDpca');
end
savefile2 = ['/Volumes/nemosine/prfsomato_fits/' currentModel sep '_myRFs_big.mat'];
save(savefile2,'myRF');

% savefile2 = ['/Volumes/nemosine/prfsomato_fits/' currentModel sep '_myBestTC.mat'];
% save(savefile2,'myBestTC');








