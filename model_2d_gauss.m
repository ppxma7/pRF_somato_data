function [thisrf, xy] = model_2d_gauss(params, stims, nSteps,dosmooth)
% model_2d_gauss - 2d gaussian model
%
%   MA's implementation of the 2d model - just breaking out for plotting
%
% ds 2022-01

%optional input args
if ieNotDefined('stims'), stims = 4; end
if ieNotDefined('nSteps'), nSteps = 3; end
if ieNotDefined('dosmooth'), dosmooth = 0; end




% calculate deltaX from nSteps
deltaX = (stims-1)/nSteps;

[stimsX, stimsY] = meshgrid(1:stims, 1:deltaX:stims);
%[stimsX, stimsY] = meshgrid(1:stims, 1:stims);

% loop through each voxel and pick rows 1..4 as x,y, stdx, stdy
thisrf = nan([size(stimsX), size(params,2)]);
for ii = 1:size(params,2)
    
    % unpack the params - refactor this?
    
    px = params(1,ii);
    py = params(2,ii);
    pstdx = params(3,ii);
    pstdy = params(4,ii);
    
    thisrf(:,:,ii) = exp(-(((stimsX-px).^2)/(2*(pstdx^2))+((stimsY-py).^2)/(2*(pstdy^2))));
    
end

if dosmooth
    % try a smoother version
    smoothFac = 20;
    %[stimsX, stimsY] = meshgrid(1:stims, 1:stims);
    
    
    inrf = thisrf;
    thisrf = nan([smoothFac, smoothFac,size(params,2)]);
    %thisrf = nan([size(stimsX), size(params,2)]);
    
    [X,Y] = meshgrid(1:size(inrf,2), 1:size(inrf,1));
    for jj = 1:size(inrf,3)
        tempRF = inrf(:,:,jj);
        F = scatteredInterpolant(X(:), Y(:), tempRF(:), 'linear');
        [U,V] = meshgrid(linspace(1,size(inrf,2),smoothFac), linspace(1,size(inrf,1),smoothFac));
        thisrf(:,:,jj) = F(U,V);
    end
    
end



if nargout == 2
    xy.stims = stims;
    xy.nSteps = nSteps;
    xy.deltaX = deltaX;
    xy.stimsX = stimsX;
    xy.stimsY = stimsY;
end

end