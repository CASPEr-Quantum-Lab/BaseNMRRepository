classdef LorentzianQuadratureFitter
    %% LorentzianQuadratureFitter class template
    %%% Brief description: %%%
    % Template for constructing an object and performing necessary
    % validation checks before fitting Lorentzian quadratures. Should speed
    % up user-time performance and improve quality of life
    %%% Author: Andrew J. Winter, Boston University %%%
    %%% Code originally written: 4/29/2025 %%%
    %%% Inputs: %%%
    % fitBoundsArray: Struct containing .startPoint, .lowerBounds, .upperBounds
    %%% Optional inputs: %%%
    % numberOfMultistarts:  Total number of multistart processes user
    %       wishes to run when performing nonlinear fit (default = 100)
    % multistartObj: Multistart object (created when calling MultiStart()
    %       function) sets multistart tolerances
    % optimizationProblemStruct: Optimization structure (contains
    %       non-linear optimization method, function tolerances, etc.)
    %%% Outputs: %%%
    % - A LorentzianQuadratureFitter object with fitting properties that
    %   able to be adjusted live (saving precious time) and input
    %   validation
    % --- Class contains many useful functions that the user can take
    %     advantage of, none of which are forcibly run (e.g. user can
    %     generate some different plots if they wish, but none of them will
    %     be forcibly run saving computation time)
    %%% Revisions: %%% 
    % 
    %%% Notes: %%%
    % %%% MATLAB dependencies: %%%
    %%% Code written with MATLAB 2022a %%%
    % REQUIRED MATLAB TOOLBOXES: 
    % - Global Optimization (for seting up the non-linear optimization)
    % SUGGESTED MATLAB TOOLBOXES:
    % - Parallel Computing (for parfor loops, speeds up code significantly)

    properties (GetAccess = public)
        NumberOfMultistarts = []; % This will be defaulted to 100 in the constructor unless otherwise provided by the user
    end
    
    properties (GetAccess = public, SetAccess = private)
        FitBoundsStruct = []; % Struct containing .startPoint, .lowerBounds, .upperBounds
        MultistartObject = []; % Forcibly set in the constructor
        OptimizationProblemStructure = []; % Forcibly set in the constructor
    end

    methods
        %%% Class constructor %%%
        function obj = LorentzianQuadratureFitter(fitBoundsStruct, numberOfMultistarts, multistartObj, optimizationProblemStruct)
            %NMRFitter Construct an instance of this class
            %   Detailed explanation goes here
            
            numberOfInputArguments = nargin; % How many input variables did the user provide?
            
            %%% Input validation
            % Check for proper number of input arguments
            if numberOfInputArguments == 0
                error("Not enough inputs provided to LorentzianQuadratureFitter: Need at least Lorentzian fit parameter start point, upper, and lower bounds.");
            end

            % Perform input validation on the fitBoundsStruct variable before assigning as a class property
            obj.verifyFitBoundsStruct(fitBoundsStruct);
            
            obj.FitBoundsStruct = fitBoundsStruct; % fitBoundsStruct is verified, assign as a class property

            if numberOfInputArguments > 1
                obj.NumberOfMultistarts = numberOfMultistarts;
            elseif numberOfInputArguments > 2
                obj.MultistartObject = multistartObj;
            elseif numberOfInputArguments > 3
                obj.OptimizationProblemStructure = optimizationProblemStruct;
            elseif numberOfInputArguments > 4
                error("Too many inputs provided to LorentzianQuadratureFitter.")
            end
            
            % Assign defaults for unprovided inputs
            if isempty(obj.NumberOfMultistarts)
                obj.NumberOfMultistarts = 100; % Default to 100
            end
            
            if isempty(obj.MultistartObject)
                obj = obj.setMultistartDefaults();
            end
            
            if isempty(obj.OptimizationProblemStructure)
                obj = obj.setOptimizationDefaults(fitBoundsStruct);
            end
            
        end
        %%% End of constructor %%%

        %%% Regular public methods
        function obj = setNumberOfMultistarts(obj, numberOfMultistarts)
            %SETNUMBEROFMULTISTARTS Summary of this method goes here
            %   Detailed explanation goes here
            obj.NumberOfMultistarts = numberOfMultistarts;
        end

        function obj = setFitBounds(obj, newFitBoundsStruct)
            %SETFITBOUNDS Verify the new fit bounds struct passes validation checks
            obj.verifyFitBoundsStruct(newFitBoundsStruct);

            obj.FitBoundsStruct = newFitBoundsStruct;
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private methods (like property-specific input validation) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function verifyFitBoundsStruct(obj, fitBoundsStruct)
            %VERIFYFITBOUNDSSTRUCT Check to make sure fitBoundsStruct is a structure not a matrix etc.
            if string(class(fitBoundsStruct)) ~= "struct"
                error("Provided fitbounds is not in the form of a struct.")
            end
            
            % Check for proper struct fields
            if ~isfield(fitBoundsStruct, 'startPoint')
                error("Provided fitBoundsStruct is missing the .startPoint field.")
            end
            if ~isfield(fitBoundsStruct, 'lowerBounds')
                error("Provided fitBoundsStruct is missing the .startPoint field.")
            end
            if ~isfield(fitBoundsStruct, 'upperBounds')
                error("Provided fitBoundsStruct is missing the .startPoint field.")
            end
            
            % Check that struct field lengths are identical
            startPointFieldLength = length(fitBoundsStruct.startPoint);
            lowerBoundsFieldLength = length(fitBoundsStruct.lowerBounds);
            upperBoundsFieldLength = length(fitBoundsStruct.upperBounds);
            tempLogicalVector = startPointFieldLength.*ones(3,1) == [startPointFieldLength; lowerBoundsFieldLength; upperBoundsFieldLength];
            if sum(tempLogicalVector) ~= 3 % There should be three "logically true" values: one for each of the fit vectors
                error("Vector lengths between the fields of the in the fitBoundsStruct are not equal.")
            end
        end

        function obj = setMultistartDefaults(obj)
            %SETMULTISTARTDEFAULTS Default multistart protocol for the nonlinear fitting (setting box tolerances etc.)
            multistartObj = MultiStart('PlotFcns',@gsplotbestf,'Display','final');
            multistartObj.TolX = 1e-3;
            multistartObj.TolFun = 1e-5;
            
            % Check if user has parallel toolbox installed before defaulting to parallel processing
            userHasParallelToolbox = ver('parallel');
            userHasParallelToolbox = exist('userHasParallelToolbox', 'var');
            if userHasParallelToolbox && obj.NumberOfMultistarts > 200
                multistartObj.UseParallel = true;
            end
            obj.MultistartObject = multistartObj;
        end

        function obj = setOptimizationDefaults(obj, FitBoundsStruct)
            %SETOPTIMIZATIONDEFAULTS Default minimization optimization problem (well-suited for Lorentzian fitting)
            optimizationProblemStruct = createOptimProblem('fmincon','x0',obj.FitBoundsStruct.startPoint, ...
                'objective',@objective,'lb',obj.FitBoundsStruct.lowerBounds,'ub',obj.FitBoundsStruct.upperBounds);
            optimizationProblemStruct.options.MaxIter = 5e3;
            optimizationProblemStruct.options.MaxFunEvals = 2e4;
            optimizationProblemStruct.options.TolX = 1e-3;
            optimizationProblemStruct.options.TolFun = 1e-6;
            optimizationProblemStruct.options.FinDiffRelStep = 1e-5;
            obj.OptimizationProblemStructure = optimizationProblemStruct;
        end
            
    end
end