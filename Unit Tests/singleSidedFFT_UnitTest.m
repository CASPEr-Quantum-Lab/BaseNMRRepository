% Simple unit test of singleSidedFFT.m function

clearvars;
close all;

debugPlots = true; % Display time and frequency domain plots
peakLocationTolerance = 1e-10; % (Hz) Tolerance to pass unit test

Fs = 2500;          % (Hz) Sampling frequency                    
T = 1/Fs;           % (s) Sampling period       
L = 1500;           % Length of signal
t = ((0:L-1)*T)';   % (s) Time vector

signalFrequencies = [50; 120]; % (Hz)
sinewaveSignal = 1*(0.7*sin(2*pi*signalFrequencies(1)*t) + sin(2*pi*signalFrequencies(2)*t)); % (mV)

sinewaveSignal = sinewaveSignal + 2*randn(size(t)); % (mV) Add pseudorandom white noise

[frequencyAxis, frequencyData] = singleSidedFFT(Fs*1e-6,sinewaveSignal);

if debugPlots
    f1 = figure;
    movegui("northwest")
    plot(t*1e3,sinewaveSignal)
    title("Signal Corrupted with Zero-Mean Random Noise")
    xlabel("time (milliseconds)")
    ylabel("mV")
end

if debugPlots
    f2 = figure;
    movegui("north")
    plot(frequencyAxis.*1e6, abs(frequencyData))
    title("Signal Corrupted with Zero-Mean Random Noise")
    xlabel("frequency (Hz)")
    ylabel("uV/sqrt(Hz)")
end

% Find the two frequency peaks
[~, peakLocations] = findpeaks(abs(frequencyData), frequencyAxis.*1e6, 'SortStr', 'descend', 'NPeaks', 2);

% Sort the initial frequencies array (this is here in case someone decides
% to change the frequencies for testing purposes)
signalFrequencies = sort(signalFrequencies, 'descend');
peakDifferences = signalFrequencies - peakLocations; % (Hz) Subtract the ordered expected frequencies from the FFT frequencies

% If any peakDifference is larger than the tolerance, assertion fails
assert(~any(peakDifferences > peakLocationTolerance), 'Peak tolerance was NOT satisfied! Max peakDifference: %0.7f', max(peakDifferences)) 