classdef irPhys < ir 
% irPhys stores general properties of the retinal patch and stores the
% rgcPhys mosaic objects in its mosaic property field. irPhys and rgcPhys
% store mosaics with RGC parameters set by fits to data from the
% Chichilnisky Lab.
% 
%       obj = irPhys(inputObj, params); 
%
% An irPhys object takes as input a bipolar object or an outerSegment
% object. The irPhys (inner retina) object stores basic properties about
% the inner retina such as the position of the simulated retinal patch.
%
% See Pillow, Jonathan W., et al. "Spatio-temporal correlations and visual
% signalling in a complete neuronal population." Nature 454.7207 (2008) and
% Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries in ON
% and OFF ganglion cells of primate retina." The Journal of Neuroscience
% 22.7 (2002). irPhys generates a mosaic of RGCs with parameters from GLM
% fits from an experiment in the Chichilnisky lab.
%
% Properties:
%     name:      animal, ir; example: 'macaque ir'
%     row:       N Stimulus row samples
%     col:       N Stimulus col samples
%     spacing:   Stimulus input spacing (um)
%     timing:    Stimulus input time step (sec)
%     eyeSide:   Left or right eye
%     eyeRadius: Position of patch in radius
%     eyeAngle:  Angle (degrees)
%     temporalEquivEcc: calculated from retinal position, see retinalLocationToTEE
%     numberTrials: number of trials for spike generation
%     mosaic: a cell array, where each cell is an rgcMosaic object,
%               which is a subclass of the ir object.
%
% Methods: set, get, compute, plot
%   see individual m-files for details.
% 
% Example:
%   innerRetinaSU = irPhys(bp, params);
% 
% See also ir.m, t_rgcCascade.m
% 
% 9/2015 JRG
% 7/2016 JRG updated

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = private, GetAccess = public)
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = irPhys(inputObj, varargin)
            % Initialize the parent class
            obj = obj@ir(inputObj, varargin{:});
            
            % Initialize ourselves by building rgcPhys mosaic objects
            obj.mosaic{1} = rgcPhys(obj, varargin{:});
            
        end
        
        % set function, see superclass method in @ir for details
        function obj = irSet(obj, varargin)
            irSet@ir(obj,varargin{:});
        end
        
        % get function, see superclass method in @ir for details
        function val = irGet(obj, varargin)
           % val = irGet(obj, varargin{:});
           val = irGet@ir(obj,varargin{:});
        end      
        
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, outersegment, varargin)
            obj = irCompute(obj, outersegment, varargin{:});
        end
        function irPlot(obj, varargin)
            irPlot@ir(obj, varargin{:});
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
