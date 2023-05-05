function [thisrf, B, xy] = model_2d_gauss(params, stims, nSteps,dosmooth)
% model_2d_gauss - 2d gaussian model
%
%   MA's implementation of the 2d model - just breaking out for plotting
%
% ds 2022-01

%optional input args
if ieNotDefined('stims'), stims = 4; end
if ieNotDefined('nSteps'), nSteps = 3; end
if ieNotDefined('dosmooth'), dosmooth = 0; end

if stims == 5
    % need to account for TW 2D Gaussian, which is actually just a 1D
    % Gaussian
    %keyboard
    stimsX = 1:stims;
    stimsY = ones(1,stims);

else
    deltaX = (stims-1)/nSteps;
    [stimsX, stimsY] = meshgrid(1:stims, 1:deltaX:stims);

end

%% loop through each voxel and pick rows 1..4 as x,y, stdx, stdy
thisrf = nan([size(stimsX), size(params,2)]);
for ii = 1:size(params,2)
    
    % unpack the params - refactor this?
    
    px = params(1,ii);
    py = params(2,ii);
    pstdx = params(3,ii);
    pstdy = params(4,ii);
    
    thisrf(:,:,ii) = exp(-(((stimsX-px).^2)/(2*(pstdx^2))+((stimsY-py).^2)/(2*(pstdy^2))));
    
end


%% SMOOTH
%smoothFac = 20;
inrf = thisrf;
if dosmooth
    % try a smoother version

    if stims == 5
        smoothFac = 20;
        % need to use interp1 instead of scattered interpolant 
        thisrf = nan([1, smoothFac,size(params,2)]);

        [xq,~] = meshgrid(linspace(1,stims,smoothFac),1);
        for jj = 1:size(inrf,3)
            tempRF = inrf(:,:,jj);
            vq1 = interp1(stimsX, tempRF ,xq);
            thisrf(:,:,jj) = vq1;
            clear vq1 tempRF
        end
        %figure, plot(stimsX,thisrf(:,:,1),'o',xq,vq1,':.')
        %figure('Position',[100 100 500 80]), imagesc(vq1)


    else
        smoothFac = 20;
        %[stimsX, stimsY] = meshgrid(1:stims, 1:stims);


        %inrf = thisrf;
        thisrf = nan([smoothFac, smoothFac,size(params,2)]);
        %thisrf = nan([size(stimsX), size(params,2)]);

        [X,Y] = meshgrid(1:size(inrf,2), 1:size(inrf,1));
        for jj = 1:size(inrf,3)
            tempRF = inrf(:,:,jj);
            F = scatteredInterpolant(X(:), Y(:), tempRF(:), 'linear');
            [U,V] = meshgrid(linspace(1,size(tempRF,2),smoothFac), linspace(1,size(tempRF,1),smoothFac));
            thisrf(:,:,jj) = F(U,V);
        end
    end

end

%% try extrapolation to fill out missing parts of RF
%interpRF = thisrf(:,:,1);
if stims ==5
    nanGrid = nan(1, 60);
else
    nanGrid = nan(60,60);
end
VV = nanGrid;
%VV(21:40,21:40) = interpRF;

% https://uk.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
% method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
for kk = 1:size(thisrf,3)
    interpRF = thisrf(:,:,kk);
    if stims == 5
        VV(1,21:40) = interpRF;
        %VV(1,51:70)= interpRF;
        %B(:,:,kk) = interp1(1:60,VV,1:60,'pchip');
        B(:,:,kk) = inpaint_nans(VV,1);
    else
        VV(21:40,21:40) = interpRF;
        B(:,:,kk) = inpaint_nans(VV,2);
    end
    %B(:,:,kk) = inpaint_nans(VV,2); % VV is matrix with holes, 2 is method
    VV = nanGrid;
end

% [X,Y] = meshgrid(1:size(B,2), 1:size(B,1));
% for jj = 1:size(inrf,3)
%     tempRF = B(:,:,jj);
%     F = scatteredInterpolant(X(:), Y(:), tempRF(:), 'linear');
%     [U,V] = meshgrid(linspace(1,size(B,2),smoothFac), linspace(1,size(B,1),smoothFac));
%     BBdwnsmpl(:,:,jj) = F(U,V);
% end


%%


if nargout == 2
    xy.stims = stims;
    xy.nSteps = nSteps;
    if stims ~= 5
        xy.deltaX = deltaX;
    end
    xy.stimsX = stimsX;
    xy.stimsY = stimsY;
end

end