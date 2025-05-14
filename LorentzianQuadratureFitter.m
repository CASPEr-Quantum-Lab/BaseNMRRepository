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
    %%% Properties: %%% (excluding inputs listed above)
    % FrequencyAxis_MHz: The frequency axis (in MHz) for the most-recent
    %       frequency-domain fit that was performed
    % ASD_uVperRtHz: The amplitude spectral density of the FID (uV/sqrt(Hz))
    %%% Public methods: %%% (excluding methods that set inputs/properties listed above)
    % fitLorentzian(): Fits the provided frequency domain data to a
    %       six-parameter Lorentzian - if real and imaginary offsets are zero,
    %       these can be forced to zero in the provided fit parameters
    %%% Revisions: %%% 
    % 
    %%% Notes: %%%
    % 
    % %%% MATLAB dependencies: %%%
    %%% Code written with MATLAB 2022b %%%
    % REQUIRED MATLAB TOOLBOXES: 
    % - Global Optimization (for seting up the non-linear optimization)
    % SUGGESTED MATLAB TOOLBOXES:
    % - Parallel Computing (for parfor loops, speeds up code significantly)

    properties (Access = public)
        NumberOfMultistarts = []; % This will be defaulted to 100 in the constructor unless otherwise provided by the user
        plotBandwidth_MHz = []; 
        fitBandwidth_MHz = [];
    end
    
    properties (GetAccess = public, SetAccess = private)
        FitBoundsStruct = []; % Struct containing .startPoint, .lowerBounds, .upperBounds
        MultistartObject = []; % Forcibly set in the constructor
        OptimizationProblemStructure = []; % Forcibly set in the constructor
        OptimizationFunction = []; % User-changeable optimization function (Lorentzian, double Lorentzian, Gaussian, etc.)
        FourierWindowDuration = [] % (us) Full duration of the time-domain window in which the Fourier transform will be taken
        FrequencyAxis_MHz = []; % (MHz) Most-recently-used frequency axis with this  fit object (SHOULD STAY SET-ACCESS PRIVATE)
        ASDmatrix_uVperRtHz = []; % (uV/sqrt(Hz)) Most-recently-used ASD with this fit object (SHOULD STAY SET-ACCESS PRIVATE)
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
                obj = obj.setOptimizationSettings();
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
        
        %%% ANDREW NOTE: THIS IS GENERALIZED! JUST NEEDS EXTRA ARGUMENT FOR
        %%% FIT TYPE
        function obj = fitLorentzian(obj, fourierWindowDuration, frequencyAxis_MHz, ASD_uVperRtHz)
            %FITLORENTZIAN Take frequency data and fit to a single
            %              Lorentzian with the fit properties of this class
            
            obj = obj.updateObjectPropertiesForFitting(fourierWindowDuration, frequencyAxis_MHz, ASD_uVperRtHz);
                        %obj.NumberOfMultistarts = 1; % DEBUGGING
            [bestFitParameters,residuals,~,~,sols] = run(obj.MultistartObject, obj.OptimizationProblemStructure, obj.NumberOfMultistarts);
            % Can we estimate error with Hessian matrix?
        end

        function obj = updateObjectPropertiesForFitting(obj, fourierWindowDuration, frequencyAxis_MHz, ASD_uVperRtHz)
            if nargin < 4
                if isempty(obj.FrequencyAxis_MHz) || isempty(obj.ASDmatrix_uVperRtHz) || isempty(obj.FourierWindowLength)
                    error("LorentzianQuadratureFitter.fitLorentzian() function supplied insufficient number of inputs!")
                end
            else
                obj.FourierWindowDuration = fourierWindowDuration;
                obj.FrequencyAxis_MHz = frequencyAxis_MHz;
                obj.ASDmatrix_uVperRtHz = ASD_uVperRtHz;
            end

            obj = obj.setOptimizationSettings();
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private methods (like property-specific input validation or necessary background functions) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        function obj = setOptimizationSettings(obj, optimizationFunction)
            %setOptimizationSettings Default minimization optimization problem (well-suited for Lorentzian fitting)
            optimizationProblemStruct = createOptimProblem('fmincon','x0',obj.FitBoundsStruct.startPoint, ...
                'objective',@objective,'lb',obj.FitBoundsStruct.lowerBounds,'ub',obj.FitBoundsStruct.upperBounds);
            optimizationProblemStruct.options.MaxIter = 5e3;
            optimizationProblemStruct.options.MaxFunEvals = 2e4;
            optimizationProblemStruct.options.TolX = 1e-3;
            optimizationProblemStruct.options.TolFun = 1e-6;
            optimizationProblemStruct.options.FinDiffRelStep = 1e-5;
            obj.OptimizationProblemStructure = optimizationProblemStruct;
            if nargin < 2
                obj.OptimizationFunction = @LorentzianQuadratureFitter.LorentzianOptimizationFunction;
            else
                obj.OptimizationFunction = optimizationFunction;
            end
            function fun = objective(fitParams)
                fprintf("test")
                %amp = p(1); gamma = p(2); f0 = p(3); phi = p(4); offset_real = p(5); offset_imag = p(6);
                
                %L_abs = (amp/sqrt(2*tD))*gamma/(2*pi)^2./((gamma/(2*pi)).^2+((f0-ff)).^2); % Real part Lorentzian
                %L_disp = (amp/sqrt(2*tD))*(ff-f0)/(2*pi)./((gamma/(2*pi)).^2+((f0-ff)).^2); % Imag part of Lorentzian
                %{
                fun = optimizationFunction(obj.FrequencyAxis_MHz, obj.ASDmatrix_uVperRtHz, ...
                                            obj.FourierWindowDuration, fitParams);
                %}
                
                fun = LorentzianQuadratureFitter.LorentzianOptimizationFunction(obj.FrequencyAxis_MHz, ...
                                            obj.ASDmatrix_uVperRtHz, obj.FourierWindowDuration, fitParams);
                %}
                %1/length(ff)*double( sum(( real(Xdata) - ( cos(phi)*L_abs + sin(phi)*L_disp + offset_real ) ).^2 + ...
                %    ( imag(Xdata) - ( sin(phi)*L_abs - cos(phi)*L_disp + offset_imag ) ).^2) ); % This is minimized
            end
        end

        function fun = objectiveFunction(fitParams)
            fprintf("test")
            %amp = p(1); gamma = p(2); f0 = p(3); phi = p(4); offset_real = p(5); offset_imag = p(6);
            
            %L_abs = (amp/sqrt(2*tD))*gamma/(2*pi)^2./((gamma/(2*pi)).^2+((f0-ff)).^2); % Real part Lorentzian
            %L_disp = (amp/sqrt(2*tD))*(ff-f0)/(2*pi)./((gamma/(2*pi)).^2+((f0-ff)).^2); % Imag part of Lorentzian
            
            fun = LorentzianQuadratureFitter.LorentzianOptimizationFunction(obj.FrequencyAxis_MHz, obj.ASD_uVperRtHz, obj.FourierWindowLength, fitParams);
            %1/length(ff)*double( sum(( real(Xdata) - ( cos(phi)*L_abs + sin(phi)*L_disp + offset_real ) ).^2 + ...
            %    ( imag(Xdata) - ( sin(phi)*L_abs - cos(phi)*L_disp + offset_imag ) ).^2) ); % This is minimized
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FITTER TYPES (e.g. Lorentzian, double Lorentzian, Gaussian, etc.) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static)
        
        function optimizationParameter = LorentzianOptimizationFunction(frequencyAxis, dataASD, fourierWindowLength, fitParameters)
            dataLength = length(dataASD);
            amplitude = fitParameters(1); gammaLinewidth = fitParameters(2);
            centerFrequency = fitParameters(3); phase = fitParameters(4);
            realOffset = fitParameters(5); imagOffset = fitParameters(6);

            realLorentzian = (amplitude/sqrt(2*fourierWindowLength))*gammaLinewidth/(2*pi)^2./((gammaLinewidth/(2*pi)).^2 ...
                            +((centerFrequency-frequencyAxis)).^2); % Real part Lorentzian
            imagLorentzian = (amplitude/sqrt(2*fourierWindowLength))*(frequencyAxis-centerFrequency)/(2*pi)./((gammaLinewidth/(2*pi)).^2 ...
                            + ((centerFrequency-frequencyAxis)).^2); % Imag part of Lorentzian
            
            optimizationParameter = 1/length(dataLength)*double( sum( ... % This is the value that is sought to be minimized
                ( real(dataASD) - ( cos(phase)*realLorentzian + sin(phase)*imagLorentzian + realOffset ) ).^2 + ...
                ( imag(dataASD) - ( sin(phase)*realLorentzian - cos(phase)*imagLorentzian + imagOffset ) ).^2) ...
                                                                    );
        end

    end
end