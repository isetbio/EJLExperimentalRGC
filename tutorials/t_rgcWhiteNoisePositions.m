% Build the mosaic with the positions from Lauren's mosaic, but give it
% default RGC sRF and tCenter and tonicDrive.

% load('/Users/james/Documents/MATLAB/isetbio misc/prosthesis_glm_fits/on_parasol_positions.mat');



%% Initialize
clear;
ieInit;

%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6;

% Stimulus length
nSteps = 5*24;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 


%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 40;
params.col = 80;
params.fov = fov;
% % params.vfov = 0.7;

%%% Grating subunit stimulus

iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;
%% Show raw stimulus for osIdentity
% % figure;
% % for frame1 = 1:size(whiteNoise.sceneRGB,3)
% %     imagesc(squeeze(whiteNoise.sceneRGB(:,:,frame1,:)));
% %     colormap gray; drawnow;
% % end
% % % close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('displayrgb');

sceneSize = sceneGet(whiteNoise.scene,'size');
retinalPatchWidth = sensorGet(whiteNoise.sensor,'width','m');
% retinalPatchHeight = sensorGet(whiteNoise.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

% timeStep = sensorGet(whiteNoise.sensor,'time interval','sec');
timeStep = (1/125)/1;
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');

%% Build RGC array

clear paramsIR innerRetina
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 4;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

model   = 'LNP';    % Computational model
innerRetina = irCreate(os,paramsIR);
% innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);

% innerRetina.mosaic{1}.mosaicSet('numberTrials',1);
% irPlot(innerRetina,'mosaic');

% % figure;

innerRetina = irSet(innerRetina,'numberTrials',1);

% filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/May25_on2/WNstim_response_OffParasol_RGC_may16.mat'];

% filenameRGC = ['C:\Users\James\Documents\GitHub\May30_onBig1\WNstim_response_OnParasol_RGC.mat'];
% save(filenameRGC, 'innerRetina');
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/WNstim_response_OnParasol_RGC.mat')
% load('/Users/james/documents/MATLAB/isetbio misc/optimal linear decoder/WNstim_response_OnParasol_RGC_may10.mat');

for blockNum = 2%:200

% clear psthNorm spikesout spikesoutM spikesoutsm whiteNoiseSmall whiteNoise iStim absorptions innerRetina

blockNum
%%% Grating subunit stimulus
% load(filenameRGC, 'innerRetina');
iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;

% filename1 = ['C:\Users\James\Documents\GitHub\May30_offBig1\WNstim_response_OffParasol_block_' num2str(blockNum) '.mat'];
% load(filename1); clear spikesoutsm;
% whiteNoise.sceneRGB = double(whiteNoiseSmall);

os = osSet(os, 'rgbData', whiteNoise.sceneRGB);

innerRetina = irCompute(innerRetina,os);

% irPlot(innerRetina, 'linear');
% irPlot(innerRetina, 'psth');

%% Look at covariance matrix
% load('/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/ws_pixiumWhiteNoise_May4.mat')
spikesout = mosaicGet(innerRetina.mosaic{1},'spikes');

% psth1 = psthstruct.psth;
% spikesout = psthstruct.spikes;

% szCells = size(psth1);
% cellCtr = 0;
% szEnd = length(psth1{1,1});
% for i2 = 1:szCells(1)
%     for j2 = 1:szCells(2)
%         cellCtr = cellCtr+1;
%         minend = min([szEnd length(psth1{i2,j2})]);
%         psthM(cellCtr,1:minend) = psth1{i2,j2}(1:minend);
%     end
% end
% 
% % rf = zeros(96,96);
% % for fr = 1:1197
% %     psthMb = (sum(psthM(17,(fr-1)*100+1:fr*100)));
% % %     rf = rf+(sum(psth1{1,2}((fr-1)*100+1:fr*100)))*(iStim.sceneRGB(:,:,fr,1));
% % %     rf = rf+(1*(iStim.sceneRGB(:,:,fr,1)));
% % end
% % figure; imagesc(rf)
% 
% for i2 = 1:szCells(1)*szCells(2)
% psthM(i2,:) = psthM(i2,:) - mean(psthM(i2,:));
% end
% 
% 
% 
% psthNorm=psthM'*diag(1./sqrt(sum(psthM'.*psthM')));
% % psthNorm = psthNorm - mean(psthNorm(:));
% psthCov = psthNorm'*psthNorm;
% % figure; imagesc(psthCov)
% % xlabel('Cell number'); ylabel('Cell number');
% % title('Covariance matrix');

% % % % % % % 
% clear spikesoutB spikesoutM
% corrfr = 1;
% for fr = corrfr+1:corrfr:1997
%     spikesoutB(:,fr-1) = sum(spikesout(:,(fr-corrfr)*100+1:fr*100));
% end
% 
% spikesoutB = spikesout;
% rf = zeros(96,96);
% for fr = 1:1997
% %     psthMb = (sum(psthM(17,(fr-1)*100+1:fr*100)));
%     rf = rf+sum(spikesout(43,(fr-1)*100+1:fr*100))*(iStim.sceneRGB(:,:,fr+1,1));
% %     rf = rf+(1*(iStim.sceneRGB(:,:,fr,1)));
% end
% figure; imagesc(rf)
% %%%%%
% % load('/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/WNstim_response_OnParasol_spikes.mat')
% for i2 = 1:szCells(1)*szCells(2)
%     spikesoutM(i2,:) = spikesoutB(i2,:) - mean(spikesoutB(i2,:));
% end
% 
% psthNorm=spikesoutM'*diag(1./sqrt(sum(spikesoutM'.*spikesoutM')));
% % psthNorm = psthNorm - mean(psthNorm(:));
% psthCov = psthNorm'*psthNorm;
% figure; imagesc(psthCov)
% % xlabel('Cell number'); ylabel('Cell number'); 
% % title('Covariance matrix');

%%%%%

% PSTHs are here:
% psthM;

% Stimulus here:
% whiteNoise.sceneRGB;

whiteNoiseSmall = uint8(squeeze(whiteNoise.sceneRGB(:,:,:,1)));
% responseSpikes = mosaicGet(innerRetina.mosaic{1},'responseSpikes');   

spikesoutsm = uint8(spikesout);
% filename1 = ['/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/May25_on2/WNstim_response_OnParasol_block_may25_' num2str(blockNum) '.mat'];
filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/WNresponse/prosthesis_on_parasol/RGC_artficial_' num2str(blockNum) '.mat'];
  
% filename1 = ['C:\Users\James\Documents\GitHub\May30_onBig1\WNstim_response_OnParasol_block_' num2str(blockNum) '.mat'];
% save(filename1, 'whiteNoiseSmall','spikesoutsm');
% toc
close
end
