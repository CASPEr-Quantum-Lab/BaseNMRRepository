% Andrew's magical test / debug script (NOT a unit test... yet)

% timestampedDirectories = ["19961008_0200_1", "19961008_0200_1"];
% dataPath = ".\testData\FIDData\";
timestampedDirectories = ["20220512_1928", "20220512_1932"];
dataPath = "D:\Github\AJWBranch\BaseNMRRepository\testData\FIDData\";
channelNumber = 1;
params.binSize = 1e3;

metaVariablesToLoad = {"cardInfo", "Nsweep"};
% Load metadata from first timestamped directory
load(strcat(dataPath, timestampedDirectories(1), "\metaFile.mat"), "-mat", metaVariablesToLoad{:});

for i = 1:1%params.NumberOfBins
    RunningAvg = 0;
    k = 1;
    for pulseIdx = (params.binSize*(i-1) + 1) : (params.binSize*i)
        SingleAvg = loadFIDDataFiles(channelNumber, dataPath, timestampedDirectories, pulseIdx);
        RunningAvg = (RunningAvg*(k-1) + SingleAvg)/k;
        k = k + 1;
        if (mod(k, 100) == 0)
            percentDone = k/params.binSize*i*100; % For displaying cmd window message loading progress
            fprintf("Loading %2.1f %%... \n", percentDone);
        end
    end
end


