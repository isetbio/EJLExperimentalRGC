function obj = irCompute(obj, inputObj, varargin)
% Computes the rgcPhys mosaic responses to an input
%
%    ir = irCompute(ir, input, varargin)
%
% Inputs:
%   ir: inner retina object
%   input - There are two types of possible input objects.
%
%      'osDisplayRGB' - frame buffer values for a spatiotemporal stimulus
%           stored in an outer segment object.
%      'bipolar' - the bipolar cell object with a signal that has gone
%           through temporal filtering and possibly spatial subunit
%           nonlinearities.
%
% Computes the linear and spike responses for all the mosaics within the
% inner retina object.  
% 
% For each mosaic, a space-time separable linear response is computed. This
% stage of the computation is stored in 'responseLinear'.  This is managed
% in irComputeLinearSTSeparable.  There is no noise added in the linear
% part.
%
% The spikes are computed irComputeSpikes routine. The spiking can have a
% random element.  So, we may run the conversion from linear to spikes
% multiple times, effectively producing spike rasters.
%
% Outputs:
%  ir: the inner retina object with responses attached to each mosaic
%
% Example:
%   ir.compute(bp);
%   irCompute(ir, bp);
%
% See also: rgcMosaic, irComputeSpikes, irComputeLinearSTSeparable
%
% JRG (c) isetbio team
% 7/2016 JRG updated

% The superclass rgcCompute carries out convolution of the linear STRF:
obj = irCompute@ir(obj, inputObj, varargin{:});

fprintf('     \n');
fprintf('Spike Generation:\n');
tic;
for cellTypeInd = 1:length(obj.mosaic)
    
    % Call Pillow code to compute spiking outputs for N trials    
    [spikeResponseFull, spikeDrive] = computeSpikesPhysLab(obj.mosaic{cellTypeInd,1});
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'responseSpikes', spikeResponseFull);
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'responseVoltage', spikeDrive);
    
    clear spikeResponseFull spikeDrive psthResponse raster psth
end
toc;




