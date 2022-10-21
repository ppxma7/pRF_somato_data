function[thisrf] = model_uncon(params, stims, nSteps, dosmooth)
% make the 16 pRF model
%stims = 4;
% stimsY = meshgrid(1:0.03:stims,1:0.03:stims);
% stimsX = transpose(stimsY);
%deltaX = (stims-1)/nSteps;
%[stimsX, stimsY] = meshgrid(1:stims, 1:deltaX:stims);



if ieNotDefined('stims'), stims = 4; end
if ieNotDefined('nSteps'), nSteps = 3; end
if ieNotDefined('dosmooth'), dosmooth = 0; end


tempGrid = params(1:16,:);
for ii = 1:size(tempGrid,2)
    thisrf(:,:,ii) = reshape(tempGrid(:,ii), [stims stims]);
    [~, thisMOI(ii)] = centerOfMass(thisrf(:,:,ii));
    %[thisfit(ii).com(jj,:), thisfit(ii).momentOfInertia(jj) ] = centerOfMass(bigmodel(:,:,jj));
end

%keyboard
% should we run the MOI now?


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



