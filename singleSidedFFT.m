%% singleSidedFFT function
%%% Brief description: %%%
% Uses sampling frequency and time-domain data to FFT and produce
% frequency domain data, uses Parseval's theorem as a validity check.
%%% Author: Andrew J. Winter, Boston University %%%
%%% Code originally written: 4/24/2025 %%%
%%% Inputs: %%%
% - samplingFrequency (MHz; DAC sampling frequency)
% - voltageDataMatrix (mV; column-wise time domain voltage data from DAC; matrices accepted)
%%% Outputs: %%%
% - frequencyAxis (MHz; corresponding frequency axis)
% - singleSidedASD (uV/sqrt(Hz); amplitude spectral density)
%%% Revisions: (GITHUB SHOULD TAKE CARE OF THIS BUT MAYBE GOOD PRACTICE TO KEEP TRACK HERE TOO) %%% 
% 
%%% Notes: %%%
% - voltageDataMatrix can be multiple sets of data, as long as its stored column-wise. This will
%   perform faster than for-looping an FFT
%%% MATLAB dependencies: %%%
%%% Code written with MATLAB 2022b %%%

function [frequencyAxis, singleSidedASD] = singleSidedFFT(samplingFrequency, voltageDataMatrix)
%% Obtain single-sided Fourier data
numberOfDataPoints = size(voltageDataMatrix, 1); % Total number of datapoints assuming data is stored column-wise
if numberOfDataPoints == 1 % If user did row-wise and has only one dataset
    numberOfDataPoints = size(voltageDataMatrix, 2); % Take the length of the row instead
end

singleSidedASD = fft(voltageDataMatrix, [], 1); % Obtain column-wise Fourier transform, becomes single-sided in the next line
singleSidedASD = singleSidedASD(1:numberOfDataPoints/2+1,:); % Include DC up to Nyquist frequency, omit the rest (this is single-sided)

% Double PSD (aka sqrt(2)*ASD) of every frequency except DC and Nyquist (since they appear twice in FFT):
singleSidedASD = [singleSidedASD(1,:); sqrt(2)*singleSidedASD(2:end-1,:); singleSidedASD(end,:)]; % (uV)
discreteFourierNormalization = sqrt(1/(samplingFrequency*numberOfDataPoints)); % (1/sqrt(Hz))
singleSidedASD = discreteFourierNormalization.*singleSidedASD; % (uV/sqrt(Hz)) To clarify units: mV/sqrt(MHz) = uV/sqrt(Hz)

%% Validity check (Parseval's theorem)
parsevalTolerance = 1e-7; % (uV^2) This is the set tolerance to pass the Parseval theorem check
powerFromFreqDomain = double(sum(abs(singleSidedASD).^2.*(samplingFrequency/numberOfDataPoints))); % (uV^2)
powerFromTimeDomain = double(mean(abs(voltageDataMatrix).^2,1)); % (uV^2)
powerMismatch = abs(powerFromTimeDomain - powerFromFreqDomain); % (uV^2)
if any(powerMismatch > parsevalTolerance) % (uV^2) Compare power mismatch with set tolerance
    error("Parseval's theorem was not obeyed: Time - Freq domain power mismatch.");
end
clear voltageDataMatrix; % voltageData is no longer useful for the rest of this function's scope: free up some memory

%% Create frequency axis
frequencyAxis = ((0:(numberOfDataPoints/2)).*(samplingFrequency/numberOfDataPoints))'; % (MHz)

end




















