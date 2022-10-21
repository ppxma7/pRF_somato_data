%% grabData 
% grab the data from the OneDrive folder
%
% 

warning('!! using an absolute path for now.. fix this ')
mypath='/Volumes/nemosine/prfsomato_fits';


% check for other users... needs fix.
if ~isfolder(mypath)
    if isunix()
        startDir = '~';
    else
        startDir = 'C:\';
    end
    mypath = uigetdir(startDir, 'Pick prfsomato_fits folder (local/onedrive)');
end

% range of data
% SUBJECTS
mysubs = {'prf1','prf2','prf3','prf4','prf6','prf7','prf8','prf9','prf10','prf11'};

% DIGITS / PD locations
digits = 2:5;
baa = 2:5;

%% set up data structure so we can loop

datafilenames = {'prf_2ddouble_dd_nov', ...
    'prf_1d_nodiag2_nov', ...
    'prf_1dt_nodiag2_nov', ...
    'prf_sixteen_nodiag2_com_nov'};

nModels = numel(datafilenames);

data = {}; % empty struct

%% load in data from summaries

for iSub = 1:length(mysubs)
    
    currentSub = mysubs{iSub};
    
    tic()
    
    % for each data file
    for iModel = 1:nModels
        currentModel = datafilenames{iModel};
        
        % basetips...
        currentFilename = [mypath filesep() currentSub filesep() 'basetips' filesep() currentModel ];
        
        data{iSub}.model(iModel).data = load(currentFilename);
        data{iSub}.model(iModel).modelname = currentModel;
        
    end
    t = toc();
    fprintf('took: %.2fs to read sub %d\n', t, iSub);
    
    % can we do BAA / ROI reading here as well?
    for iDigit = 1:numel(digits)
        currentDigit = digits(iDigit);
        
        for iBaa = 1:numel(baa)
            currentBaa = baa(iBaa);
            % just pick one model (ROI indeces are all the same
            % load([mypath mysubs{iSub} '/baa/' 'D' num2str(iDig) 'pd' num2str(iBaa) '_16.mat'],'prf_overlays');
            currentROIFilename = [mypath filesep() currentSub filesep() '/baa/' sprintf('D%dpd%d_16.mat', currentDigit, currentBaa)];
            data{iSub}.roiIndex(iDigit, iBaa) = load(currentROIFilename, 'idx');
        end
    end
    
end

