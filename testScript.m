% Andrew's magical test / debug script (NOT a unit test... yet)
addpath(".\CASPErHelpers\")
clearvars;
close all;

debugPlots = false;

% Fit parameter arrays for single Lorentzian go as:
% [time domain amplitude (mV) , linewidth gamma (MHz) , 
%  center frequency (MHz) , phase (rad) ,
%  real offset (uV/sqrt(Hz)) , imag offset (uV/sqrt(Hz) )
% ]

chosenSweepIdxs = [2];
fourierWindowStart = 7; % (us) Start with respect to end of pulse - so this can go negative for excitation spectra!
fourierWindowDuration = 1.5e3; % (us) Duration, starting at fourierWindowStart, to Fourier transform

initialGuess = [0.038, 0.03, 4.6, pi, 0, 0];

gammaDeviation = 0.9; % Low and high bound deviation of gamma (multiplied for lower bound, division for upper bound)
fitBoundsStruct = struct( ...
            "startPoint", initialGuess, ...
            "lowerBounds", [-1e4*initialGuess(1), initialGuess(2)*0.9, initialGuess(3)*0.95, 0, -10, 0], ...
            "upperBounds", [1e4*initialGuess(1), initialGuess(2)/0.9, initialGuess(3)*1.05, 2*pi, 10, 0] ...
                        );
fitterObject = LorentzianQuadratureFitter(fitBoundsStruct);

%timestampedDirectories = ["19961008_0200_1", "19961008_0200_2"];
%dataPath = ".\testData\FIDData\";
% totalNumberOfMeasurementsToLoad: Total number of datafiles to analyze PER directory
% (e.g. user selects 5, and dir1 has 5 files dir2 has 10 files, only first 5 of
% each are loaded)
%totalNumberOfMeasurementsToLoad = 5; % ONE "measurement" = ONE .bin file
%dataBinSize = 5; % Desired number of data files to bin-average together
timestampedDirectories = ["20220512_1928"];%, "20220512_1932"];
totalNumberOfMeasurementsToLoad = 500; % ONE "measurement" = ONE .bin file
dataPath = "C:\Users\Andrew\Desktop\MATLAB Scripts\Quantum Lab\207PbGen2Data\";
dataBinSize = 500;
channelNumber = 1;

metaVariablesToLoad = {"cardInfo", "Nsweep", "signal"};
% Load metadata from first timestamped directory
load(strcat(dataPath, timestampedDirectories(1), "\metaFile.mat"), "-mat", metaVariablesToLoad{:});
fivePercentProgressBarNotifier = int16(dataBinSize*0.05); % This is used for a progress bar when loading
if fivePercentProgressBarNotifier == 0
    fivePercentProgressBarNotifier = 1;
end
numberOfDataBins = totalNumberOfMeasurementsToLoad/dataBinSize; % Calculate the total number of data bins expected
if numberOfDataBins < 1 % User input error
    error("Number of bins is less than one, check the 'totalNumberOfMeasurementsToLoad' and 'dataBinSize' variables!")
end
% Take data and average in bins of size dataBinSize
for binIdx = 1:numberOfDataBins
    timeDomainAvg = 0;
    k = 1;
    for pulseIdx = (dataBinSize*(binIdx-1) + 1) : (dataBinSize*binIdx)
        singleAvg = loadFIDDataFiles(channelNumber, dataPath, timestampedDirectories, pulseIdx);
        % Running average of binned datafiles
        timeDomainAvg = (timeDomainAvg*(k-1) + singleAvg)/k; % (mV)
        k = k + 1;
        if (mod(k-1, fivePercentProgressBarNotifier) == 0)
            percentDone = (k-1)/dataBinSize*binIdx*100; % For displaying cmd window message loading progress
            fprintf("Loading / binning %2.0f%%... \n", percentDone);
        end
    end
end
clear singleAvg; % Garbage clean-up
sampleRate_MHz = cardInfo.setSamplerate.*1e-6; % (Hz --> MHz)

fullTimeAxis = (0:1:(size(timeDomainAvg,1)-1))'./cardInfo.setSamplerate.*1e6; % (s --> us)

%% Janos shit

pulseDuration = signal{1}.Duration(1)*1e6; % pulse duration [us]
pulseStart = signal{1}.Start(1)*1e6; % pulse start [us]
pulseStart = 1; %Alex

if isempty(chosenSweepIdxs) 
    chosenSweepIdxs = 1:Nsweep;
end
try
    pulseRepetitionTime = signal{1}.Start(2) - signal{1}.Start(1);
catch
    pulseRepetitionTime = 0;
end

% end of pulse
params.pulseEnd = pulseDuration + pulseStart; % Alex

% set window
fourierWindowStart = pulseDuration + pulseStart + fourierWindowStart; %[us] Fourier window starts

%%% Intermediate Andrew time
[minValueStart, startTimeIdx] = min(abs(pulseStart - fullTimeAxis));
[minValueEnd, endTimeIdx] = min(abs(pulseStart + fourierWindowDuration - fullTimeAxis));
newTimeAxis = fullTimeAxis(startTimeIdx:endTimeIdx);
%%% End intermediate Andrew time

% FFT time <-> index conversion
numberOfPointsUpToFourierStart = round(sampleRate_MHz*fourierWindowStart);
numberOfPointsFromTriggerToStart = numberOfPointsUpToFourierStart + round(sampleRate_MHz*fourierWindowDuration); %nF + 2^(nextpow2(round(Fs * tD))-1);
windowLengthFromTrigger = numberOfPointsFromTriggerToStart/sampleRate_MHz; % NOT USED
totalNumberOfTimeIdxs = numberOfPointsFromTriggerToStart - numberOfPointsUpToFourierStart;
tTot = totalNumberOfTimeIdxs / sampleRate_MHz;
windowedTimeAxis = ((numberOfPointsUpToFourierStart + 1 : numberOfPointsFromTriggerToStart)/sampleRate_MHz); %time in s
%timeDomainWindow = windowedTimeAxis - (pulseDuration + pulseStart); % SEEMINGLY NOT USED
newFrequencyAxis = sampleRate_MHz*(0:(totalNumberOfTimeIdxs/2))/totalNumberOfTimeIdxs; % (MHz)

%% End Janos shit

% Convert time domain data to single-sided frequency domain
[frequencyAxis, ASDmatrix] = singleSidedFFT(sampleRate_MHz, timeDomainAvg(startTimeIdx:endTimeIdx,:));

if debugPlots % Plot time domain data
    f1 = figure;
    movegui('northwest')
    plot(fullTimeAxis, timeDomainAvg)
    xline(windowedTimeAxis(1), '--k', 'LineWidth', 1.5)
    xline(windowedTimeAxis(end), '--k', 'LineWidth', 1.5)
    xlabel('time (\mus)')
    ylabel('voltage at digitizer input (mV)')
    CASPErFigSettings(f1);
    
    %[frequencyAxis, ASDmatrix] = singleSidedFFT(sampleRate_MHz, timeDomainAvg);
    f2 = figure;
    movegui('north')
    plot(frequencyAxis, abs(ASDmatrix), '.')
    xlabel('frequency (MHz)')
    ylabel('voltage at digitizer input (uV/\surd(Hz))')
    CASPErFigSettings(f2);
end

%% Load the data into the quadrature fitter and fit
chosenFrequencyDatasets = ASDmatrix(:, chosenSweepIdxs);

fitterObject = fitterObject.fitLorentzian(fourierWindowDuration, frequencyAxis, chosenFrequencyDatasets);





