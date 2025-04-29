%% PREVIOUSLY KNOWN AS LoadDataFiles.m 
%% Details:
%% This code loads FID datafiles (collected as single precision data from FID MODE in the GUI, NOT AXION MODE which uses int16)
%% v1.0: baseline version. Currently have to manually switch between *_avg.bin and regular .bin FID data types

%% Need:
% channelNumber (scalar)
% timestampedDirectories (vector of directories) <-- these contain all the pulse measurements
% dataPath (character array) <-- path to the directory containing data "timestampedDirectories" directories
% measurementNumber (scalar) <-- specific data file idx, particularly
%               useful in binning (i.e. bins of 50 will properly count as
%               vector of 0-49, 50-99, etc.)
% 

function [averagedData, numAvgVector] = loadFIDDataFiles(channelNumber, dataPath, timestampedDirectories, pulseNumber, useAveragedData)
    
    if nargin == 4 % If user does not specify whether or not to use averaged data
        useAveragedData = true; % Default to using averaged data (if exists)
    end

    metaVariablesToLoad = {'cardInfo.AI', 'segmentSize', 'Nsweep'};
    % Load metadata from first timestamped directory
    load(strcat(dataPath, timestampedDirectories(1), "\metaFile.mat"), '-mat', metaVariablesToLoad);
    
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
        load(strcat(fullpath, "\metaFile.mat"), '-mat', 'numAvg');
        
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
