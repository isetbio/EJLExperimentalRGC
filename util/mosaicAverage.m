function mosaicAverageGLM = mosaicAverage(mosaicGLM)
% Computes the average properties of a mosaicGLM fitted to data from a
% Chichilnisky Lab experiment. 
% 
%       mosaicAverageGLM = mosaicAverage(mosaicGLM); 
%             [called only from @rgcPhys/mosaicLoadAverage.m]
% 
% The fit parameters from each cell in the mosaic GLM are pulled out and
% averaged. These parameters include the spatial receptive field, the
% difference of Gaussians fit to the spatial receptive field, the temporal
% impulse response, the tonic drive and the nonlinearity.
% 
% See also @rgcPhys/mosaicLoadAverage.m
% 
% 7/2016 JRG (c) isetbio team

% Get linear filters from each cell
for i = 1:length(mosaicGLM)
    
    % Get spatial filter
    if isfield(mosaicGLM{i}.linearfilters.Stimulus,'space_rk1')
        sRFtemp = mosaicGLM{i}.linearfilters.Stimulus.space_rk1;
        sRF(i,:,:) = sRFtemp;
        
        % Measure max, min and mean values
        maxsRF(i)  = max(sRFtemp(:));
        minsRF(i)  = min(sRFtemp(:));
        meansRF(i) = mean(sRFtemp(:));
    end
    
    % Get temporal filter
    if isfield(mosaicGLM{i}.linearfilters.Stimulus,'time_rk1')
        
        tCtemp = mosaicGLM{i}.linearfilters.Stimulus.time_rk1;
        
        % Measure max, min and mean values
        tC(i,:) = tCtemp;
        maxtC(i) = max(tCtemp);
        mintC(i) = min(tCtemp);
        meantC(i) = mean(tCtemp);
    end
    
    % Get tonic drive filter
    if isfield(mosaicGLM{i}.linearfilters.Stimulus,'tonicDrive')            
        tonicD(i,:) = mosaicGLM{i}.linearfilters.Stimulus.tonicDrive{i,1}(:);
    elseif isfield(mosaicGLM{i}.linearfilters,'TonicDrive')
        tonicD(i,:) = mosaicGLM{i}.linearfilters.TonicDrive.Filter;
    else
        tonicD(i,:) = 0;
    end
    
    % Get nonlinearity
    if isfield(mosaicGLM{i},'model');
        nlcoeffs(i,:) = mosaicGLM{i}.model.Coefficients.Estimate;
    end
    
    % Get DoG fit
    if isfield(mosaicGLM{i},'stafit')
        if isfield(mosaicGLM{i}.stafit,'center_sd_x')
            sd_x(i,:) = mosaicGLM{i}.stafit.center_sd_x;
            sd_y(i,:) = mosaicGLM{i}.stafit.center_sd_y;
        end
    end
    
end

%%
% % Plot spatial RF
% meanRF = mean(sRF);
% sRFrs = reshape(meanRF,size(sRFtemp,1),size(sRFtemp,2));
% figure; surf(sRFrs);

% Average spatial RFs
% Center each RF on max value before averaging
oldSize = size(squeeze(sRF(i,:,:)),1);
newSize = size(sRFtemp,1)*3+1;
meanRF = zeros(newSize);
meanCtr = zeros(newSize);

for i = 1:length(mosaicGLM)
    % Find max value
    [maxPR(i) maxPRind(i)] = max(max(abs(squeeze(sRF(i,:,:))),[],1));
    [maxPC(i) maxPCind(i)] = max(max(abs(squeeze(sRF(i,:,:))),[],2));
   
    % Center max value
    xv = [(floor(newSize/2)+1) - floor(oldSize/2) : (floor(newSize/2)+1) + floor(oldSize/2)] - (maxPCind(i) - floor(oldSize/2));
    yv = [(floor(newSize/2)+1) - floor(oldSize/2) : (floor(newSize/2)+1) + floor(oldSize/2)] - (maxPRind(i) - floor(oldSize/2));  

    % Check sign
    sRFi = squeeze(sRF(i,:,:));
    [m1,m2] = max(abs(sRFi(:)));  
    signmult = sign(sRFi(m2));    
    tC(i,:) = signmult*tC(i,:);
    
    % Add to mean
    meanRF(xv,yv) = meanRF(xv,yv) + signmult*squeeze(sRF(i,:,:));
    meanCtr(xv,yv) = meanCtr(xv,yv) + ones(size(squeeze(sRF(i,:,:))));
end

% % % Plot mean spatial RF
% figure;
% for i = 1:length(mosaicGLM)
%     sRFi = squeeze(sRF(i,:,:));
%     [m1,m2] = max(abs(sRFi(:)));  
%     signmult = sign(sRFi(m2));
%     subplot(121);
%     surf(squeeze(signmult*sRF(i,:,:)));
%     subplot(122);
%     plot(tC(i,:));
% end
% figure; surf(meanRF./length(mosaicGLM)); shading flat

%% Store average values for isetbio format
cv = [(floor(newSize/2)+1) - floor(oldSize/2) : (floor(newSize/2)+1) + floor(oldSize/2)] - 1;
mosaicAverageGLM.linearfilters.Stimulus.space_rk1 = meanRF(cv,cv)./length(mosaicGLM); % meanAvg(cv,cv);
mosaicAverageGLM.linearfilters.Stimulus.time_rk1 = mean(tC);
mosaicAverageGLM.linearfilters.Stimulus.tonicDrive = mean(tonicD);

if isfield(mosaicGLM{i},'model');
    mosaicAverageGLM.modelavg = mean(nlcoeffs);

    % Need to change this
    mosaicAverageGLM.model = mosaicGLM{1}.model;
end

if isfield(mosaicGLM{i},'stafit')
    mosaicAverageGLM.sd = [sqrt(mean(sd_x(find((sd_x~=0)&(sd_x<6)&(sd_x>0.4))).^2)) sqrt(mean(sd_y(find((sd_y~=0)&(sd_y<6)&(sd_y>0.4))).^2))];
end
