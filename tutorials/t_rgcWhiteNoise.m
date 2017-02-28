% t_rgcSelectiveRecon
% 
% Simulate white noise resopnse
% 
% Simulate stimulus reconstruction

%% Initialize 
clear
% ieInit;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 3;%:2    % Choose On Parasol (1) or Off Parasol (2)
stimulusTestI = 2;%:2     % Choose WN test stimulus (1) or NSEM test stimulus (2)
    
% Switch on the conditions indices
% Experimental dataset
switch experimentI
    case 1; experimentID = '2013-08-19-6';
    otherwise; error('Data not yet available');
end
% The other experimental data will be added to the RDT in the future.

% Stimulus: white noise or natural scene movie with eye movements
switch stimulusTestI
    case 1; stimulusTest = 'WN';
    case 2; stimulusTest = 'NSEM';
end

% Cell type: ON or OFF Parasol
switch cellTypeI
    case 1; cellType = 'prosthesis selective';
    case 2; cellType = 'prosthesis off parasol';
    case 3; cellType = 'prosthesis on parasol';
end

%% Load stimulus movie and fit/spiking data using RemoteDataToolbox

% Loads the appropriate movie and spiking data for the experimental
% conditions.
% [testmovie, xval_mosaic] =  loadDataRGCFigure2(experimentI,stimulusTestI,cellTypeI);


% Length of WN movie is 1200, take nFrames to limit natural movie to same length
% nFrames = 3600; 
% testmovieshort = double(testmovie.matrix(:,:,1:nFrames)); 

%%
for blockNum = 1:2%:20
 
clear spikesout whiteNoiseSmall whiteNoise iStim innerRetina innerRetinaSpikes
 
blockNum

%% Set horizontal field of view
fov = 8;
 
% Stimulus length
nSteps = 24000*5;
 nFrames = nSteps;
% Activation curve
 
% Spatial activation of electrode
 
% Electrode PWM 
  
% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
paramsWN.nSteps = nSteps;
paramsWN.row = 40;
paramsWN.col = 80;
paramsWN.fov = fov;
% % paramsWN.vfov = 0.7;
 
 
iStim = ieStimulusBinaryWhiteNoise(paramsWN);
absorptions = iStim.sensor;
whiteNoise = iStim;
testmovieshort = iStim.sceneRGB;

% clear whiteNoiseSmall
% filename0 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/WNresponse/prosthesis_off_parasol/WNstim_response_ProsOFFParasol_block_' num2str(blockNum) '.mat'];
% load(filename0);
% testmovieshort = whiteNoiseSmall;
%% Upsample movie stimulus for os computation    
frRS = 1;
testmovieRS = zeros(size(testmovieshort,1),size(testmovieshort,2),frRS*size(testmovieshort,3));
for frnum = 1:nFrames
    for frrep = 1:frRS
        testmovieRS(:,:,(frnum-1)*frRS+frrep) = testmovieshort(:,:,frnum);
    end
end

%% Generate outer segment object for GLM from RGB scene data

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', 1+double(testmovieshort));

%% Generate RGC object for simulated GLM prediction of response
% Set the parameters for the inner retina RGC mosaic. For the inner retina
% type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
% effect, because those are dependent on the properties of the retinal
% piece used in the Chichilnisky Lab experiment.

% Set parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = 12; 
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning to allow looping
params.experimentID = experimentID; % Experimental dataset
params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol

% params.cellIndices = 10;%:10%:118;

% Create object
innerRetina = irPhys(os1, params);
nTrials = 1;
innerRetina.mosaic{1} = innerRetina.mosaic{1}.set('numberTrials',nTrials);
%% Compute the inner retina response

innerRetina = irCompute(innerRetina, os1);

% Get the PSTH from the object
% innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');

innerRetinaSpikes = mosaicGet(innerRetina.mosaic{1},'spikes');


whiteNoiseSmall = uint8(squeeze(iStim.sceneRGB(:,:,:,1)));
% responseSpikes = mosaicGet(innerRetina.mosaic{1},'responseSpikes');   
 imagesc(whiteNoiseSmall(:,:,20)); drawnow
mosSize = innerRetina.mosaic{1}.get('mosaicsize');
spikesoutsm = zeros(mosSize(1),24000);
spikesoutsm(:,1:size(innerRetinaSpikes,3)) = uint8(squeeze(innerRetinaSpikes));
% filename1 = ['/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/May25_on2/WNstim_response_OnParasol_block_may25_' num2str(blockNum) '.mat'];
% filename1 = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_block_' num2str(blockNum) '.mat'];
filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/WNresponse/prosthesis_on_parasol/WNstim_response_ProsONParasol_stixds4_short6_' num2str(blockNum) '.mat'];
save(filename1, 'whiteNoiseSmall','spikesoutsm');

end