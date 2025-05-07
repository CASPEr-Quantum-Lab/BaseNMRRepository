%% Feedback Fitter Script
% Andrew Mattson, 4/08/2025
% This script loads in the (FFT) results of backgroundAverager and fits the data to the noise equations as defined in the Feedback Bible.
% First, with no feedback to extract parameters that are unchanged by feedback. Then, with feedback to get the rest.

%%
clearvars;
close all;


%% Load in data and variables

% change to location of SPN2024_HardcodedParameters.mat file
load('/Users/andrewmattson/Desktop/Research/SPN2024_HardcodedParameters.mat',...
     'Ltp','Ltp_err','L1_0','L1_0_err','L2','L2_err','M12_0','M12_0_err','C_0','C_0_err',...
     'Lsq','Lsq_err','SQUIDs0076_Min','SQUIDs0076_Mfb','Rf','GBW','lowFreqSQUIDtransfer',...
     'extFluxResistance','calToneAmplitude','calToneAttenuation','oldExpectedCalToneRMS','warmGainFactorA0','warmGainFactorA0_err','frequencyCorrectionFactorb','frequencyCorrectionFactorb_err','bA0','bA0_err',...
     'coilRadius','centralArmL','turnNumber',...
     'H1Tneel','H1gyromag','F19Tneel','F19gyromag',...
     'kB','hbar','mu0', ...
     "datasets", "df",...
     "meanFitF0","meanFitF0_err","meanFitOffsets","meanFitOffsets_err","meanFitCalVals","meanFitCalVals_err",...
     "cap","cap_err","Lc","Lc_err","Lp","L1","M12","M12_err","G0","G0_err"...
     )
SQUIDcurrentCalibration = lowFreqSQUIDtransfer*bA0; % Overall gain (SQUID + warm gain + frequency correction) as described in Supplement Section 10.2.  (V/pA)*(V/V).  Use digitizerVoltagePSD/SQUIDcurrentCalibration^2 to get S_{II}


%% does the data contain signal tones for SNR?
SNR = 0; % 0 means no, 1 means yes

%% Load in variables

% first path should be data with no feedback, the rest are the ones with different feedback settings 
% !! Change path base to match file locations !!
dataPathBase = '/Users/andrewmattson/Desktop/Research/FeedbackTesting/20250314_';
dataPaths = {'1922', '2024', '1918', '1357'};  % timestamps produced by background averager
feedback_gains = {'0', '15.8 V/V', '31.6 V/V', '? V/V (1357)'}; % for plot labeling
% dataPathBase = '/Users/andrewmattson/Desktop/Research/FeedbackTesting/20250416_';
% dataPaths = {'1742', '1739'};
% feedback_gains = {'0', '34.8 V/V'};

fitRange = [3.68,3.72]*1e6; % Hz
plotRange = [3.68,3.72]*1e6; % Hz

for i = 1:length(dataPaths)
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

    % if no SNR
    if SNR == 0
        fitPSD{i} = temp(i1:i2)/SQUIDcurrentCalibration^2; %pA^2/Hz
        plotPSD{i} = temp(i1p:i2p)/SQUIDcurrentCalibration^2; %pA^2/Hz
    end

    % if SNR is desired
    if SNR == 1
        % find signal index
        sig_freq = 3.766e6; 
        [~, fit_sig_idx] = min(abs(fitFreq - sig_freq));
        [~, plot_sig_idx] = min(abs(plotFreq - sig_freq));
        fit_sig_idx = fit_sig_idx;
        plot_sig_idx = plot_sig_idx;

        % get PSDs
        fitPSDtemp(:,i) = temp(i1:i2)/SQUIDcurrentCalibration^2; %pA^2/Hz
        plotPSDtemp(:,i) = temp(i1p:i2p)/SQUIDcurrentCalibration^2; %pA^2/Hz

        % remove it from the data so it doesn't skew fit
        fit_noise_idx = (1:length(fitPSDtemp(:,i))) ~= fit_sig_idx;
        plot_noise_idx = (1:length(plotPSDtemp(:,i))) ~= plot_sig_idx;

        % make new arrays
        fitFreq_signal = fitFreq(fit_sig_idx);
        fitFreq = fitFreq(fit_noise_idx);
        plotFreq_signal = plotFreq(plot_sig_idx);
        plotFreq = plotFreq(plot_noise_idx);
        fitPSD_signal{i} = fitPSDtemp(fit_sig_idx,i); %pA^2/Hz
        fitPSD{i} = fitPSDtemp(fit_noise_idx,i); %pA^2/Hz
        plotPSD_signal{i} = plotPSDtemp(plot_sig_idx,i); 
        plotPSD{i} = plotPSDtemp(plot_noise_idx,i);

        plotPSDfull{i} = plotPSDtemp(:,i);
    end
end

%% Do the fit with no feedback

% Fit parameters
stptnfb = [1e4,3.7e6,218.3,0.05];
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
nfb_cf1_ = fit(fitFreq,fitPSD{1},nfb_ft1_,nfbfo1_);
nfb_cf1_


%% Do the fit with feedback
if length(dataPaths) > 1
    for i = 2:length(dataPaths)
        disp('...... Fitting')   
        
        % Fit parameters:
        stpt = [1e4,nfb_cf1_.b,nfb_cf1_.g, 0.5, pi,   nfb_cf1_.Simp, 0.007];
        lb =   [-1e10,nfb_cf1_.b,nfb_cf1_.g, 0,    0,    nfb_cf1_.Simp, 0];
        ub =   [1e10,nfb_cf1_.b,nfb_cf1_.g, 100000,  2*pi, nfb_cf1_.Simp, 1];
        
        fo1_ = fitoptions('method','NonlinearLeastSquares','Robust','off','MaxFunEvals',1500,'MaxIter',10000,...
        'Startpoint',stpt,'Lower',lb,'Upper',ub,'TolFun',fttol,'TolX',1e-12); %,'Display','iter'

        % Full model with feedback corrections
        ft1_ = fittype('(a*(g - 2*af*a*cos(phi))/(2*pi))/(((g - 2*a*af*cos(phi))/2)^2 + (x-b+af*a*sin(phi))^2) + Simp*((g/2)^2 + (x-b)^2)/((g - 2*a*af*cos(phi)/2)^2 + (x-b+a*af*sin(phi))^2) + Snn',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a','b','g','af','phi','Simp','Snn'});
    
        % Do fit
        cf1_{i} = fit(fitFreq,fitPSD{i},ft1_,fo1_);
        cf1_{i}

        % extract SNR (if in the data)
        if SNR == 1

        end

    end
end


%% Calculate residuals and make plots

% set up colors
col = hsv(length(dataPaths));

if SNR == 1
    % find SNRs
    snr_list = zeros(length(dataPaths));
    snr_list(1) = plotPSD_signal{1} / feval(nfb_cf1_,plotFreq_signal);
    for i = 2:length(dataPaths)
        snr_list(i) = plotPSD_signal{i} / feval(cf1_{i},plotFreq_signal);
    end
end

if SNR == 0
    residuals{1} = (plotPSD{1} - feval(nfb_cf1_,plotFreqfull))./feval(nfb_cf1_,plotFreq);
    % plot no feedback
    f1 = figure('Color',[1,1,1], 'Units','normalized', 'Position',[0.25,0.25,0.35,0.5]); box on; hold on
    plot(1e-3*plotFreq,plotPSD{1},'-','Color',col(1,:),'DisplayName','0 V/V')
    plot(1e-3*plotFreq,feval(nfb_cf1_,plotFreq),'-','LineWidth',2,'Color',col(1,:)/1.5,'HandleVisibility', 'off')
    
    % plot others
    if length(dataPaths) > 1
        for i = 2:length(dataPaths)
            residuals{i} = (plotPSD{i} - feval(cf1_{i},plotFreqfull))./feval(cf1_{i},plotFreq);
            plot(1e-3*plotFreq,plotPSD{i},'-','Color',col(i,:),'DisplayName',feedback_gains{i})
            plot(1e-3*plotFreq,feval(cf1_{i},plotFreq),'-','LineWidth',2,'Color',col(i,:)/1.5,'HandleVisibility', 'off');
        end
    end
end

if SNR == 1
    residuals{1} = (plotPSDfull{1} - feval(nfb_cf1_,plotFreqfull))./feval(nfb_cf1_,plotFreqfull);
    % plot no feedback
    f1 = figure('Color',[1,1,1], 'Units','normalized', 'Position',[0.25,0.25,0.35,0.5]); box on; hold on
    plot(1e-3*plotFreqfull,plotPSDfull{1},'-','Color',col(1,:),'DisplayName',strcat('0 V/V, SNR=', num2str(snr_list(1))))
    plot(1e-3*plotFreq,feval(nfb_cf1_,plotFreq),'-','LineWidth',2,'Color',col(1,:)/1.5,'HandleVisibility', 'off')
    plot(1e-3*plotFreq_signal, plotPSD_signal{1}, 'ko', 'MarkerSize', 4, 'MarkerFaceColor',col(1,:),'HandleVisibility', 'off');
    
    % plot others
    if length(dataPaths) > 1
        for i = 2:length(dataPaths)
            residuals{i} = (plotPSDfull{i} - feval(cf1_{i},plotFreqfull))./feval(cf1_{i},plotFreqfull);
            plot(1e-3*plotFreqfull,plotPSDfull{i},'-','Color',col(i,:),'DisplayName',strcat(feedback_gains{i},', SNR=', num2str(snr_list(i))))
            plot(1e-3*plotFreq,feval(cf1_{i},plotFreq),'-','LineWidth',2,'Color',col(i,:)/1.5,'HandleVisibility', 'off');
            plot(1e-3*plotFreq_signal, plotPSD_signal{i}, 'ko', 'MarkerSize', 4, 'MarkerFaceColor',col(i,:),'HandleVisibility', 'off');
        end
    end
end

% other plot settings
xlabel("frequency (kHz)")
ylabel("S_{II} @ SQUID Input (pA^2/Hz)")
grid on; grid minor;
xlim([-inf, inf]);
ylim([-inf, inf])
ax1 = gca;
set(gca,'FontSize',20,'Fontweight','normal','FontName','Arial')
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
xlabel("frequency (kHz)")
ylabel("Residual (data - fit)/(fit)")
grid on; grid minor;
xlim([-inf, inf]);
ylim([-inf, inf])
ax1 = gca;
set(gca,'FontSize',20,'Fontweight','normal','FontName','Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca,'YScale','log')
legend('Location','best')


%% the end