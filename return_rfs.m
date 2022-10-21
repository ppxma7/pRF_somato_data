function  [thisrf, thisr2, this_desc] = return_rfs(data, iSub, iDigit, iBaa, iModel, N, data_desc)
% return_rfs - return RFs from the summary data structure
%
% eg: 
%     iSub = 9;
%     iDigit = 1;
%     iBaa = 1;
%     iModel = 1;
% 
%     thisrf = return_rfs(data, iSub, iDigit, iBaa, iModel)
%
%     montage(thisrf)
%
% ds 2022-03-29

% top K rfs to return
if ieNotDefined('N'), N = 16; end


% data description
if ieNotDefined('data_desc')
    data_desc.mysubs = {'prf1','prf2','prf3','prf4','prf6','prf7','prf8','prf9','prf10','prf11'};
    data_desc.digits = 2:5;
    data_desc.baa = 2:5;
    data_desc.datafilenames = {'prf_2ddouble_dd_nov', ...
                        'prf_1d_nodiag2_nov', ...
                        'prf_1dt_nodiag2_nov', ...
                        'prf_sixteen_nodiag2_com_nov'};
end

% get info into variables
currentSub = data_desc.mysubs{iSub};
currentDigit = data_desc.digits(iDigit);
currentBaa = data_desc.baa(iBaa);
currentModel = data_desc.datafilenames{iModel};


% current data
currentD = data{iSub}.model(iModel).data.dd;

% ROI - subset that's in the particular ROI -- same as in linearCoords
currentROI = data{iSub}.roiIndex(iDigit, iBaa).idx;
 
% using the coords
% @TODO: notes to add to README...
%
%   data{iSub}.model(iModel).data.dd 
% .   contains linearCoords, which are the matlab-linear coords for the 
% .   rawCoords that may contain NANs
%
% .  to make sure the linear coords index into the correct place in
% vectors, data... take the nan's out of  
% params, r2, mytSeries, mymodelResp


% these index into the 3d space (and don't contain NANs)
linearCoords = currentD.linearCoords;

% there are in 3d, and may still contain NANs
invalidCoords =  isnan( currentD.rawCoords(1,:)) | ...
                       isnan( currentD.rawCoords(2,:)) | ...
                       isnan( currentD.rawCoords(3,:));
                    
theR2 = currentD.r2(~invalidCoords);

% all DATA / params
extractedParams = currentD.params(:, ~invalidCoords);

% INTERSECT to find the linearCoords (and their indices) to map into data
% [C,IA,IB] = intersect(A,B)
% C = A(IA) and C = B(IB)
% common, BIG then SMALL
[commonVoxels, iInData, iInROI] = intersect(linearCoords, currentROI );

% in ROI versions (not all subject)
roiR2 = theR2(iInData);

% topN = prctile(roiR2, 90); % a cutoff value
% topNidx = roiR2 >= topN; % an index in terms of data vec cooordinates
% matchingLinearIdx = iInData(topNidx);

[topN topNidxInRoi] = maxk(roiR2, N); % a cutoff value
% find maxk in ROI... from there in linear index and from there back down
% to the vec data indices
[~, matchingLinearIdx, ~] = intersect(linearCoords, currentROI(iInROI(topNidxInRoi)));

% pick the correct model function:

if contains(currentModel, '_2ddouble_')
    modelFunc = @model_2d_gauss;   
elseif contains(currentModel, '_1d_')
    modelFunc = @model_1d_gauss;
elseif contains(currentModel, '_1dt_')
    modelFunc = @model_1dT_gauss; % TRANSPOSE
elseif contains(currentModel, '_sixteen_')
    modelFunc = @model_uncon;
else
    error('(uhoh!) modelFunc cannot be defined - check your inputs')
end
    fprintf('using modelFunc = @%s\n', func2str(modelFunc) );
    

% get the RF shapes
thisrf = modelFunc(extractedParams(:, matchingLinearIdx));
% thisrf = permute(thisrf, [2 1 3]); % matrix to image order

thisr2 = topN;

this_desc.matchingLinearIdx = matchingLinearIdx;
this_desc.N = N; 
this_desc.currentSub = currentSub;
this_desc.currentDigit = currentDigit;
this_desc.currentBaa = currentBaa;
this_desc.currentModel = currentModel;

end