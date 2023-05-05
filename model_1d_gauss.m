function [thisrf, B] = model_1d_gauss(params, stims, nSteps,doT,dosmooth)
% model_1d_gauss - 1d gaussian model
%
%   MA's implementation of the 1d model - just breaking out for plotting
%
% ds 2022-01

%optional input args
if ieNotDefined('stims'), stims = 4; end
if ieNotDefined('nSteps'), nSteps = 3; end
if ieNotDefined('doT'), doT = 0; end
if ieNotDefined('dosmooth'), dosmooth = 0; end


% make the 1D Gaussian pRF model

X = linspace(1,stims,nSteps);

for ii = 1:size(params,2)
    pone = [1         params(1,ii) params(5,ii) 0];
    ptwo = [params(9,ii) params(2,ii) params(6,ii) 0];
    pthr = [params(10,ii) params(3,ii) params(7,ii) 0];
    pFour = [params(11,ii) params(4,ii) params(8,ii) 0];
    
    Z = [gauss(pone,X)', gauss(ptwo,X)', gauss(pthr,X)', gauss(pFour,X)'];
    
    if doT == 1
        thisrf(:,:,ii) = transpose(Z);
    else
        thisrf(:,:,ii) = Z;
    end
    
    
end

if dosmooth
    smoothFac = 20;    
    inrf = thisrf;
    thisrf = nan([smoothFac, smoothFac,size(params,2)]);
    [X,Y] = meshgrid(1:size(inrf,2), 1:size(inrf,1));
    for jj = 1:size(inrf,3)
        tempRF = inrf(:,:,jj);
        F = scatteredInterpolant(X(:), Y(:), tempRF(:), 'linear');
        [U,V] = meshgrid(linspace(1,size(inrf,2),smoothFac), linspace(1,size(inrf,1),smoothFac));
        thisrf(:,:,jj) = F(U,V);
    end
    
end
%% extrapolation

nanGrid = nan(60, 60);
VV = nanGrid;
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
    VV(21:40,21:40) = interpRF;
    B(:,:,kk) = inpaint_nans(VV,2); % VV is matrix with holes, 2 is method
    VV = nanGrid;
end

%% return to 20x20?
[X,Y] = meshgrid(1:size(B,2), 1:size(B,1));
for jj = 1:size(inrf,3)
    tempRF = B(:,:,jj);
    F = scatteredInterpolant(X(:), Y(:), tempRF(:), 'linear');
    [U,V] = meshgrid(linspace(1,size(B,2),smoothFac), linspace(1,size(B,1),smoothFac));
    BBdwnsmpl(:,:,jj) = F(U,V);
end





end









