%% Feedback Fitter Script
%%% Brief description: %%%
% This script loads in the (FFT) results of backgroundAverager and
% fits the data to the noise equations as defined in the Feedback Bible.
% First, with no feedback to extract parameters that are unchanged
% by feedback. Then, with feedback to get the rest. 
%%% Author: Andrew H. Mattson, Boston University %%%
%%% Code originally written: 5/5/2025 %%%
%%% Inputs: %%%
% Noise data avgd_data.mat file with no feedback
% Other noise data avgd_data.mat files with feedback
%%% Outputs: %%%
% Fit parameters for noise equations
% Residuals from fit
% Plots of fits and residuals
%%% Revisions: %%%
% 5/8/2025 (commited to AMattson branch)
% - Updated fit equation
% - Updated header to match template
% - Changed paths to relative paths
% - Changed all strings to character arrays
% - Preallocated arrays before looping
% - Other code cleaning and housekeeping
%%% MATLAB dependencies: %%%
%%% Code written with MATLAB 2022a %%%
% REQUIRED MATLAB TOOLBOXES: 
% - Symbolic Math (for legendreP function)


%%
clearvars;
close all;


%% Load in data and variables

% change to location of SPN2024_HardcodedParameters.mat file
load('./SPN2024_HardcodedParameters.mat',...
     'Ltp','Ltp_err','L1_0','L1_0_err','L2','L2_err','M12_0','M12_0_err','C_0','C_0_err',...
     'Lsq','Lsq_err','SQUIDs0076_Min','SQUIDs0076_Mfb','Rf','GBW','lowFreqSQUIDtransfer',...
     'extFluxResistance','calToneAmplitude','calToneAttenuation','oldExpectedCalToneRMS','warmGainFactorA0','warmGainFactorA0_err','frequencyCorrectionFactorb','frequencyCorrectionFactorb_err','bA0','bA0_err',...
     'coilRadius','centralArmL','turnNumber',...
     'H1Tneel','H1gyromag','F19Tneel','F19gyromag',...
     'kB','hbar','mu0', ...
     'datasets', 'df',...
     'meanFitF0','meanFitF0_err','meanFitOffsets','meanFitOffsets_err','meanFitCalVals','meanFitCalVals_err',...
     'cap','cap_err','Lc','Lc_err','Lp','L1','M12','M12_err','G0','G0_err'...
     )
SQUIDcurrentCalibration = lowFreqSQUIDtransfer*bA0; % Overall gain (SQUID + warm gain + frequency correction) as described in Supplement Section 10.2.  (V/pA)*(V/V).  Use digitizerVoltagePSD/SQUIDcurrentCalibration^2 to get S_{II}


%% does the data contain signal tones for SNR calculation?
isthereSNR = 1; % 0 means no, 1 means yes

%% Load in variables

% first path should be data with no feedback, the rest are the ones with different feedback settings 
% !! Change path base to match file locations !!
dataPathBase = './FeedbackTesting/20250314_';
dataPaths = {'1922', '2024', '1918', '1357'};  % timestamps produced by background averager
feedback_gains = {'0', '15.8 V/V', '31.6 V/V', '? V/V (1357)'}; % for plot labeling
% dataPathBase = '/Users/andrewmattson/Desktop/Research/FeedbackTesting/20250416_';
% dataPaths = {'1742', '1739'};
% feedback_gains = {'0', '34.8 V/V'};

fitRange = [3.69,3.71]*1e6; % Hz
plotRange = [3.69,3.71]*1e6; % Hz

numberOfDataPaths = length(dataPaths);

%% preallocate arrays before loop
% load in example data to get size of arrays
load([dataPathBase,dataPaths{1},'/00001/avgd_data.mat'],'df','startInd','v')
Freq = df*(startInd-2)+df*[1:length(v{1})].';
[~,i1] = min(abs(Freq-fitRange(1)));
[~,i2] = min(abs(Freq-fitRange(2)));
[~,i1p] = min(abs(Freq-plotRange(1)));
[~,i2p] = min(abs(Freq-plotRange(2)));
fitFreq = Freq(i1:i2);
plotFreq = Freq(i1p:i2p);
plotFreqfull = Freq(i1p:i2p);
temp = v{1};

% set up arrays
fitPSD = cell(numberOfDataPaths,1);
plotPSD = cell(numberOfDataPaths,1);
lengthOfFit = length(temp(i1:i2)/SQUIDcurrentCalibration^2);
lengthOfPlot = length(temp(i1p:i2p)/SQUIDcurrentCalibration^2);
fitPSDtemp = zeros(lengthOfFit, numberOfDataPaths);
plotPSDtemp = zeros(lengthOfPlot, numberOfDataPaths);
fitPSDSignal = cell(numberOfDataPaths,1);
plotPSDSignal = cell(numberOfDataPaths,1);
plotPSDfull = cell(numberOfDataPaths,1);


%% loop through data files and perform data extraction and trimming
for i = 1:numberOfDataPaths
    load([dataPathBase,dataPaths{i},'/00001/avgd_data.mat'],'df','startInd','v')
    Freq = df*(startInd-2)+df*[1:length(v{1})].';
    
    % trim frequency arrays
    [~,i1] = min(abs(Freq-fitRange(1)));
    [~,i2] = min(abs(Freq-fitRange(2)));
    [~,i1p] = min(abs(Freq-plotRange(1)));
    [~,i2p] = min(abs(Freq-plotRange(2)));
    fitFreq = Freq(i1:i2);
    plotFreq = Freq(i1p:i2p);
    plotFreqfull = Freq(i1p:i2p);
    temp = v{1};

    % if no SNR calculation is wanted
    if isthereSNR == 0
        fitPSD{i} = temp(i1:i2)/SQUIDcurrentCalibration^2; %pA^2/Hz
        plotPSD{i} = temp(i1p:i2p)/SQUIDcurrentCalibration^2; %pA^2/Hz
    end

    % if SNR calculation is desired
    if isthereSNR == 1
        % find signal index
        sigFreq = 3.766e6; % Hz
        [~, fitSigIdx] = min(abs(fitFreq - sigFreq));
        [~, plotSigIdx] = min(abs(plotFreq - sigFreq));

        % get PSDs
        fitPSDtemp(:,i) = temp(i1:i2)/SQUIDcurrentCalibration^2; %pA^2/Hz
        plotPSDtemp(:,i) = temp(i1p:i2p)/SQUIDcurrentCalibration^2; %pA^2/Hz

        % remove it from the data so it doesn't skew fit
        fitNoiseIdx = (1:length(fitPSDtemp(:,i))) ~= fitSigIdx;
        plotNoiseIdx = (1:length(plotPSDtemp(:,i))) ~= plotSigIdx;

        % make new arrays
        fitFreqSignal = fitFreq(fitSigIdx);
        fitFreq = fitFreq(fitNoiseIdx);
        plotFreqSignal = plotFreq(plotSigIdx);
        plotFreq = plotFreq(plotNoiseIdx);
        fitPSDSignal{i} = fitPSDtemp(fitSigIdx,i); %pA^2/Hz
        fitPSD{i} = fitPSDtemp(fitNoiseIdx,i); %pA^2/Hz
        plotPSDSignal{i} = plotPSDtemp(plotSigIdx,i); 
        plotPSD{i} = plotPSDtemp(plotNoiseIdx,i);

        plotPSDfull{i} = plotPSDtemp(:,i);
    end
end

%% Do the fit with no feedback

% Fit parameters
stptnfb = [1e4,3.7e6,218.3,0.05]; % [amplitude (pA^2), center frequency (Hz), gamma linewidth (Hz), SQUID imprecision noise (pA^2)]
lbnfb =   [1e1,3.60e6,100,0];
ubnfb =   [1e8,3.8e6,1e4,1e3];
fttol = 1e-20;
nfbfo1_ = fitoptions('method','NonlinearLeastSquares','Robust','off','MaxFunEvals',1500,'MaxIter',10000,...
    'Startpoint',stptnfb,'Lower',lbnfb,'Upper',ubnfb,'TolFun',fttol,'TolX',1e-12); %,'Display','iter'

% Lorentzian model with integral a and FWHM g
nfb_ft1_ = fittype('(a*g/(2*pi))/((g/2)^2 + (x-b)^2) + Simp*((g/2)^2 + (x-b)^2)/((g/2)^2 + (x-b)^2)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b','g','Simp'});

% Do fit
nfb_cf1_ = fit(fitFreq,fitPSD{1},nfb_ft1_,nfbfo1_)


%% Do the fit with feedback
if numberOfDataPaths > 1
    for i = 2:numberOfDataPaths
        disp('...... Fitting')   
        
        % Fit parameters:
        stpt = [1e4, nfb_cf1_.b, nfb_cf1_.g, nfb_cf1_.b, nfb_cf1_.g, nfb_cf1_.Simp, 0.007]; % [amplitude (pA^2), center frequency (Hz), gamma linewidth (Hz), feedback adjusted gamma (Hz), feedback adjusted center freq (Hz), SQUID imprecision noise (pA^2), constant noise offset (pA^2)]
        lb =   [-1e10, nfb_cf1_.b, nfb_cf1_.g, 3e6, 0, nfb_cf1_.Simp, 0];
        ub =   [1e10, nfb_cf1_.b, nfb_cf1_.g, 4e6, 1e5, nfb_cf1_.Simp, 1];
        
        fo1_ = fitoptions('method','NonlinearLeastSquares','Robust','off','MaxFunEvals',1500,'MaxIter',10000,...
        'Startpoint',stpt,'Lower',lb,'Upper',ub,'TolFun',fttol,'TolX',1e-12); %,'Display','iter'

        % Full model with feedback corrections
        ft1_ = fittype('(a*gt/(2*pi))/((gt/2)^2 + (x-bt)^2) + Simp*((g/2)^2 + (x-b)^2)/((gt/2)^2 + (x-bt)^2) + Snn',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a','b','g','bt','gt','Simp','Snn'});
    
        % Do fit
        cf1_{i} = fit(fitFreq,fitPSD{i},ft1_,fo1_);
        cf1_{i} % for some reason it doesn't print parameters in command line unless I do it on a separate line
    end
end


%% Calculate residuals and make plots

% set up colors
col = hsv(numberOfDataPaths);

if isthereSNR == 1
    % find SNRs
    snr_list = zeros(numberOfDataPaths);
    snr_list(1) = plotPSDSignal{1} / feval(nfb_cf1_,plotFreqSignal);
    for i = 2:numberOfDataPaths
        snr_list(i) = plotPSDSignal{i} / feval(cf1_{i},plotFreqSignal);
    end
end

residuals = cell(numberOfDataPaths,1);
if isthereSNR == 0
    residuals{1} = (plotPSD{1} - feval(nfb_cf1_,plotFreqfull))./feval(nfb_cf1_,plotFreq);
    % plot no feedback
    f1 = figure('Color',[1,1,1], 'Units','normalized', 'Position',[0.25,0.25,0.35,0.5]); box on; hold on
    plot(1e-3*plotFreq,plotPSD{1},'-','Color',col(1,:),'DisplayName','0 V/V')
    plot(1e-3*plotFreq,feval(nfb_cf1_,plotFreq),'-','LineWidth',2,'Color',col(1,:)/1.5,'HandleVisibility', 'off')
    
    % plot others
    if numberOfDataPaths > 1
        for i = 2:numberOfDataPaths
            residuals{i} = (plotPSD{i} - feval(cf1_{i},plotFreqfull))./feval(cf1_{i},plotFreq);
            plot(1e-3*plotFreq,plotPSD{i},'-','Color',col(i,:),'DisplayName',feedback_gains{i})
            plot(1e-3*plotFreq,feval(cf1_{i},plotFreq),'-','LineWidth',2,'Color',col(i,:)/1.5,'HandleVisibility', 'off');
        end
    end
end

if isthereSNR == 1
    residuals{1} = (plotPSDfull{1} - feval(nfb_cf1_,plotFreqfull))./feval(nfb_cf1_,plotFreqfull);
    % plot no feedback
    f1 = figure('Color',[1,1,1], 'Units','normalized', 'Position',[0.25,0.25,0.35,0.5]); box on; hold on
    plot(1e-3*plotFreqfull,plotPSDfull{1},'-','Color',col(1,:),'DisplayName',strcat('0 V/V, SNR=', num2str(snr_list(1))))
    plot(1e-3*plotFreq,feval(nfb_cf1_,plotFreq),'-','LineWidth',2,'Color',col(1,:)/1.5,'HandleVisibility', 'off')
    plot(1e-3*plotFreqSignal, plotPSDSignal{1}, 'ko', 'MarkerSize', 4, 'MarkerFaceColor',col(1,:),'HandleVisibility', 'off');
    
    % plot others
    if numberOfDataPaths > 1
        for i = 2:numberOfDataPaths
            residuals{i} = (plotPSDfull{i} - feval(cf1_{i},plotFreqfull))./feval(cf1_{i},plotFreqfull);
            plot(1e-3*plotFreqfull,plotPSDfull{i},'-','Color',col(i,:),'DisplayName',strcat(feedback_gains{i},', SNR=', num2str(snr_list(i))))
            plot(1e-3*plotFreq,feval(cf1_{i},plotFreq),'-','LineWidth',2,'Color',col(i,:)/1.5,'HandleVisibility', 'off');
            plot(1e-3*plotFreqSignal, plotPSDSignal{i}, 'ko', 'MarkerSize', 4, 'MarkerFaceColor',col(i,:),'HandleVisibility', 'off');
        end
    end
end

% other plot settings
xlabel('frequency (kHz)')
ylabel('S_{II} @ SQUID Input (pA^2/Hz)')
grid on; grid minor;
xlim([-inf, inf]);
ylim([-inf, inf])
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'YScale','log')
legend('Location','best')
hold off;


%% plot residuals 

f2 = figure('Color',[1,1,1], 'Units','normalized', 'Position',[0.25,0.25,0.35,0.5]); box on; hold on
for i = 1:length(residuals)
    plot(1e-3*plotFreqfull,residuals{i},'-','Color',col(i,:),'DisplayName',feedback_gains{i})
end

% other plot settings
xlabel('frequency (kHz)')
ylabel('Residual (data - fit)/(fit)')
grid on; grid minor;
xlim([-inf, inf]);
ylim([-inf, inf])
set(gca,'XMinorTick','on','YMinorTick','on')
legend('Location','best')


%% Use CASPEr fig settings
CASPErFigSettings(f1);
CASPErFigSettings(f2);


%% the end