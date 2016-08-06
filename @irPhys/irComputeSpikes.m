function obj = irComputeSpikes(obj, varargin)
% Computes the spiking rgcPhys mosaic responses to an input
%
%    ir = irComputeSpikes(ir, input, varargin)
%
% Inputs:
%   ir: inner retina object
%   input - the irPhys object with the linear responses of each RGC already
%       computed.
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
%   irComputeSpikes(ir);
%
% See also: rgcMosaic, irCompute, irComputeLinearSTSeparable
%
% JRG (c) isetbio team

% The superclass rgcCompute carries out convolution of the linear STRF:
% obj = irCompute@ir(obj, outerSegment, varargin{:});

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




