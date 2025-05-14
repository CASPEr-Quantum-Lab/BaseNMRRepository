%% loadFIDDataFiles
%%% Brief description: %%%
% Function that, given a directory path, will use the metadata in a/the
%       provided directory(ies) and load + convert the raw binary data files
%       from FID experiments into mV datapoints
%%% Author: Andrew J. Winter, Boston University %%%
%%% Code originally written: 4/29/2025 %%%
%%% Inputs: %%%
% channelNumber: The base-1 channel number of the data from the digitizer
%       to be loaded (e.g. 1, 2, 3, or 4)
% dataPath: The base path leading up to, but EXCLUDING the timestamped
%       datadirectories (e.g. "F:\FIDData\")
% timestampedDirectories: A list of timestamped data directories
%       (e.g. ["20220512_1928", "20220512_1932"])
% pulseNumber: In an FID experiment with more than one measurement, this
%       designates which measurement number (i.e. .bin files are listed
%       similar to 0_000001.bin, 0_000002.bin, ...,
%       0_0000NN.bin, ..., 0_NNNNNN.bin. Ignoring the leading zeros, the
%       last number is the "pulseNumber" or measurement number.)
%       ***(e.g. 0_000193.bin corresponds to pulseNumber = 193)***
% useAveragedData: (true/false) Use averaged data if it exists
%       (It will check if avgd_data.mat exists, and if not, will go through
%       the number of binary files the user originally chose)
%%% Outputs: %%%
% - averagedData: (mV) The averaged time domain data of the FID signal
% - numAvgVector: Vector containing the number of files that were averaged,
%       in same order as the list of the timestampedDirectories vector, per
%       directory (this is useful for running average calculations)
%%% Revisions: %%% 
% 
%%% Notes: %%%
% 
% %%% MATLAB dependencies: %%%
%%% Code written with MATLAB 2022b %%%

function [averagedData, numAvgVector] = loadFIDDataFiles(channelNumber, dataPath, timestampedDirectories, pulseNumber, useAveragedData)
    
    if nargin == 4 % If user does not specify whether or not to use averaged data
        useAveragedData = true; % Default to using averaged data (if exists)
    end

    metaVariablesToLoad = {"cardInfo", "Nsweep"};
    % Load metadata from first timestamped directory
    load(strcat(dataPath, timestampedDirectories(1), "\metaFile.mat"), "-mat", metaVariablesToLoad{:});
    
    % Store the number of measurements per directory being averaged
    % together (allows for proper averaging of directories with mismatched
    % number of datafiles)
    numAvgVector = zeros(1,length(timestampedDirectories)); % Preallocate array

    % Bit to mV conversion factor
    digitizerConversion = cardInfo.AI.setRange/2^(cardInfo.AI.resolution-1).*ones(1,cardInfo.setChannels); % (mV/bit)

    for dirIdx = 1:length(timestampedDirectories)
        % Load the number of averaged files per directory/directory for
        % proper running average (in case of mis-matched directory lengths)
        fullpath = strcat(dataPath, timestampedDirectories(dirIdx));
        load(strcat(fullpath, "\metaFile.mat"), "-mat", "numAvg");
        
        numAvgVector(1,dirIdx) = numAvg; % Store this directories numAvg for later
        averageDataFilExists = exist(strcat(fullpath,'\',num2str(channelNumber-1,'%01.0f'),'_avgd','.bin'), 'file');
        if averageDataFilExists && useAveragedData
            filename = strcat(fullpath,'\',num2str(channelNumber-1,'%01.0f'),'_avgd','.bin'); % Use averaged data
        else
            filename = strcat(fullpath,'\',num2str(channelNumber-1,'%01.0f'),'_',num2str(pulseNumber,'%06.0f'),'.bin'); % Load individual data
        end
        
        % Open the data file(s) and read into matlab
        try
            fId0 = fopen(filename,'r');
            digitizerData = fread(fId0,'single'); % (bits)
            fileClose = fclose(fId0);
        catch
            error(strcat("Hey stupid, could not locate and open specified file: ", filename, ". Did you select to use averaging with unaveraged data???"))
        end
        
        % Convert binary data to mV using spectrum card bits/mV conversion
        digitizerData(:,dirIdx) = single(digitizerData) * digitizerConversion(channelNumber); % (mV)
    end
    
    % Find the absolute total number of files to be averaged between all directories
    totalNumberOfAveragedFiles = sum(numAvgVector);
    % Multiply each dataset by their corresponding number of averaged data,
    % then sum them and divide by total number of averages to produce this
    % set of multiple-directory-average
    % ANDREW WANTS TO DOUBLE CHECK THAT THIS MAKES SENSE
    averagedData = sum(digitizerData.*numAvgVector, 2)./(totalNumberOfAveragedFiles); % Take into account number of total measurements performed in the running average;
    if Nsweep > 1 % Reshape array into a matrix containing one column for each sweep frequency
        averagedData = reshape(averagedData, length(averagedData)/Nsweep, Nsweep); % Each measurement includes 1 sequence           
    end
end
