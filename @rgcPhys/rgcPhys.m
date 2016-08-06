classdef rgcPhys % < rgcMosaic
% Generates an rgcPhys mosaic object. rgcPhys is used for loading RGC
% parameters from a physiology experiment in the Chichilnisky Lab.
%
% The RGC models are detailed in Chichilnisky & Kalmar, J. Neurosci (2002);
% Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J. Neurosci (2005);
% and Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli,
% Nature (2008).
% 
% The computational model implemented here relies on code by
% <http://pillowlab.princeton.edu/code_GLM.html Pillow>, which is
% distributed under the GNU General Public License.
%
% rgcPhys is not a subclass of rgcMosaic, but is similar to rgcGLM in many
% respects. It is called when creating a new Phys model for an
% inner retina object.  Typically we get here from the inner retina object
% via a call
%
%   ir.mosaicCreate('model','phys','type','your type goes here')
% 
% See also: irPhys.m, t_rgcCascade.m
%
% Example:
%   innerRetinaSU = irPhys(bp, params); % rgcPhys called internally
%   innerRetinaSU.mosaic{1}
%
% 9/2015 JRG    
% 7/2016 JRG updated
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = protected, GetAccess = public)
        experimentID;
        cellType;           % on parasol, off parasol, on midget, off midget, sbc
        cellID;             % id number of cell from fitting procedure
        stimulusFit;        % type of stim used to fit parameters, white noise (WN) or natural scene with eye movements (NSEM)
        stimulusTest;       % type of stim used to test accuracy of model, WN or NSEM
        rfDiameter;         % the diameter of the RF in input pixels
        rfDiaMagnitude;     % 1 STD magnitude of spatial RF, for making movies of response in input pixels
        cellLocation;       % spatial center of each RF  in input pixels
        sRFcenter;          % spatial RF of the center on the receptor grid in units of conditional intensity
        sRFsurround;        % spatial RF of the surround in units of conditional intensity
        tCenter;            % temporal impulse response of the center in units of conditional intensity
        tSurround;          %    and of the surround (1 ms timing by default) in units of conditional intensity
        responseLinear;     % Store the linear response after convolution in units of conditional intensity
        
        generatorFunction;  % the mapping from linear signal to probability of spiking
        nlResponse;         % generatorFunctio(responseLinear)
        numberTrials;       % number of trials over which spikes are comptued
        spikeResponse;      % Store the spike times of the responses
        
        postSpikeFilter;    % the filter imposing inhomogeneity on Poisson spiking
        couplingFilter;     % the filter coupling the responses of nearby cells
        couplingMatrix;     % the matrix of connections between nearby neurons
        tonicDrive;         % DC term for linear response
        responseRaster;     % Store the rasters of the spike response
        responsePsth;       % Store the PSTH of the spike response
        responseVoltage;    % Store the "membrane voltage" output of the GLM
        responseSpikes;     % Store the spike times of the responses
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcPhys(rgc, varargin)
            obj = obj.initialize(rgc, varargin{:});            
        end
        
        % set function, see for details
        function obj = set(obj, varargin)
            % obj = set@rgcMosaic(obj, varargin);
            obj = mosaicSet(obj, varargin{:});
        end
        
        % get function, see for details
        function val = get(obj, varargin)
           val = mosaicGet(obj, varargin{:});
        end
        
        % find the relevant spatial region of input for each RGC
        [stimX, stimY, offset] = obj.stimPositions(xcell,ycell)
      
        % load a set of fitted RGC parameters from an experiment in the
        % Chichilnisky lab. see @rgcPhys/glmLoad.m.
        mosaicGLM = obj.glmLoad(cellType);
        
        % convert the fitted parameter structure into an rgcPhys object
        % see @rgcPhys/mosaicLoadExperimental.m
        obj = obj.mosaicLoadExperimental(mosaicGLM, cellType, varargin);
        
        % convert the fitted parameters structure into an rgcPhys object
        % with parameters extrapolated to a new (usually more foveal)
        % eccentricity.
        % see @rgcPhys/mosaicLoadAverage.m
        obj = obj.mosaicLoadAverage(mosaicGLM, cellType, varargin);
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        obj = initialize(obj, varargin);
    end
    
end
