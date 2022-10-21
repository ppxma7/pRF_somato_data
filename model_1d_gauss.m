function [thisrf] = model_1d_gauss(params, stims, nSteps,doT,dosmooth)
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



end