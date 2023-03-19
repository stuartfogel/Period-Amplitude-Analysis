function PAA_calculation_default_all()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate slow wave data summary (from .csv format output)
%
% Info:     Opens csv output from 'PAA_standalone' pipeline and computes
%           summary stats.
%
% Download: https://github.com/stuartfogel/PAA
%
% Copyright, Sleep Well. https://www.sleepwellpsg.com
%
% Date:     November 26, 2020
%
% Moficiations:
%   16-02-21: SF, fixed compatability for writetable introduced in R2019b.
%   16-02-21: SF, calculate positive and negative HW separately.
%   30-08-22: SF, minor variable name change for compatability fix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User defined parameters
PARAM.channels = {'Fz','Cz','Pz'}; % channels to extract slow wave info. Default = {'Fz','Cz','Pz','Oz'}.
PARAM.stages = {'N2','N3'}; % channels to extract slow wave info. Default = {'N2','N3'}.

%% Specify filename(s)
% you can manually specify filenames here, or leave empty for pop-up
PARAM.pathname = ''; % directory where csv file is located
PARAM.filename = ''; % names of csv file
PARAM.resultDir = ''; % directory to save results output

%% Open interface to select *.csv file(s)
if isempty(PARAM.filename)
    disp('Please select the file to process.');
    [filename,pathname] = uigetfile(   {'*.csv', 'comma-separated file (*.CSV)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Choose files to process', ...
        'Multiselect', 'on');
end

% check the filename(s)
if isequal(filename,0) % no files were selected
    disp('User selected Cancel')
    return;
else
    if ischar(filename) % only one file was selected
        filename = cellstr(filename); % put the filename in the same cell structure as multiselect
    end
end

PARAM.filename = filename;
PARAM.pathname = pathname;

%% Output directory
if isempty(PARAM.resultDir)
    disp('Please select a directory in which to save the results.');
    resultDir = uigetdir('', 'Select the directory in which to save the results');
end

PARAM.resultDir = resultDir;

clear filename pathname resultDir

disp('Processing selected files...')

%% Separate into slow waves from specified channel, stage & filename (defined above)
for nfile = 1:length(PARAM.filename)
    % read in raw data to table format
    tables{nfile} = readtable([char(PARAM.pathname) char(PARAM.filename(nfile))]);
    % take only slow waves from specified channel
    for nch = 1:length(PARAM.channels)
        perChannelidx = strcmp(PARAM.channels(nch),tables{nfile}.Channel);
        perChanData{nfile,nch} = tables{nfile}(perChannelidx,:);
        % take only slow waves from specified sleep stages
        for nstage = 1:length(PARAM.stages)
            perStageidx = strcmp(PARAM.stages(nstage),perChanData{nfile,nch}.SleepStage);
            perStageData{nfile,nch,nstage} = perChanData{nfile,nch}(perStageidx,:);
        end
    end
end

clear tables perChannelidx perStageidx nfile nch nstage

warning( 'off', 'MATLAB:xlswrite:AddSheet' );

%% Export per channel per stage per N to excel
for nch = 1:length(PARAM.channels)
    for nstage = 1:length(PARAM.stages)
        allTypeFilename = [PARAM.resultDir filesep char(PARAM.channels(nch)) '_' char(PARAM.stages(nstage)) '.xlsx']; % all slow waves
        for nfile = 1:length(PARAM.filename)
            % export to excel individual data
            sheetname = char(PARAM.filename(nfile));
            sheetname = sheetname(1:end-4);
            sheetname(isspace(sheetname)) = [];
            % write to xlsx
            if length(sheetname)>30
                error('Input file name too long. Please rename files with shorter names.')
            else
                writetable(struct2table(PARAM),allTypeFilename,'Sheet','PARAM')
                writetable(perStageData{nfile,nch,nstage},allTypeFilename,'Sheet',sheetname)
            end
        end
    end
end

clear allTypeFilename sheetname

%% Export summary data to excel

% create empty cell structures
NumberAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
DurationAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
FrequencyAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
PeakAmplitudeAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
AreaAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
AverageAmplitudeAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
UpSlopeAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
DownSlopeAllType{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];

NumberPos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
DurationPos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
FrequencyPos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
PeakAmplitudePos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
AreaPos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
AverageAmplitudePos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
UpSlopePos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
DownSlopePos{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];

NumberNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
DurationNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
FrequencyNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
PeakAmplitudeNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
AreaNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
AverageAmplitudeNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
UpSlopeNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];
DownSlopeNeg{length(PARAM.filename),length(PARAM.channels),length(PARAM.stages)} = [];

% calculate means for each channel and stage
for nfile = 1:length(PARAM.filename)
    for nch = 1:length(PARAM.channels)
        for nstage = 1:length(PARAM.stages)
            sheetname = char(PARAM.filename(nfile));
            sheetname = sheetname(1:end-4);
            ID{nfile} = sheetname;
            % all slow waves
            NumberAllType{nfile,nch,nstage} = sum(length([perStageData{nfile,nch,nstage}.Latency;NumberAllType{nfile,nch,nstage}]));
            DurationAllType{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Duration;DurationAllType{nfile,nch,nstage}]);
            FrequencyAllType{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Frequency;FrequencyAllType{nfile,nch,nstage}]);
            PeakAmplitudeAllType{nfile,nch,nstage} = mean([abs(perStageData{nfile,nch,nstage}.PeakAmp);PeakAmplitudeAllType{nfile,nch,nstage}]);
            AreaAllType{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Area;AreaAllType{nfile,nch,nstage}]);
            AverageAmplitudeAllType{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.AvgAmp;AverageAmplitudeAllType{nfile,nch,nstage}]);
            UpSlopeAllType{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.UpSlope;UpSlopeAllType{nfile,nch,nstage}]);
            DownSlopeAllType{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.DownSlope;DownSlopeAllType{nfile,nch,nstage}]);
            % postive half waves (take values with PeakAmp > 0)
            NumberPos{nfile,nch,nstage} = sum(length([perStageData{nfile,nch,nstage}.Latency(perStageData{nfile,nch,nstage}.PeakAmp>0);NumberPos{nfile,nch,nstage}]));
            DurationPos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Duration(perStageData{nfile,nch,nstage}.PeakAmp>0);DurationPos{nfile,nch,nstage}]);
            FrequencyPos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Frequency(perStageData{nfile,nch,nstage}.PeakAmp>0);FrequencyPos{nfile,nch,nstage}]);
            PeakAmplitudePos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.PeakAmp(perStageData{nfile,nch,nstage}.PeakAmp>0);PeakAmplitudePos{nfile,nch,nstage}]);
            AreaPos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Area(perStageData{nfile,nch,nstage}.PeakAmp>0);AreaPos{nfile,nch,nstage}]);
            AverageAmplitudePos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.AvgAmp(perStageData{nfile,nch,nstage}.PeakAmp>0);AverageAmplitudePos{nfile,nch,nstage}]);
            UpSlopePos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.UpSlope(perStageData{nfile,nch,nstage}.PeakAmp>0);UpSlopePos{nfile,nch,nstage}]);
            DownSlopePos{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.DownSlope(perStageData{nfile,nch,nstage}.PeakAmp>0);DownSlopePos{nfile,nch,nstage}]);
            % negative half waves (take values with PeakAmp < 0)
            NumberNeg{nfile,nch,nstage} = sum(length([perStageData{nfile,nch,nstage}.Latency(perStageData{nfile,nch,nstage}.PeakAmp<0);NumberNeg{nfile,nch,nstage}]));
            DurationNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Duration(perStageData{nfile,nch,nstage}.PeakAmp<0);DurationNeg{nfile,nch,nstage}]);
            FrequencyNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Frequency(perStageData{nfile,nch,nstage}.PeakAmp<0);FrequencyNeg{nfile,nch,nstage}]);
            PeakAmplitudeNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.PeakAmp(perStageData{nfile,nch,nstage}.PeakAmp<0);PeakAmplitudeNeg{nfile,nch,nstage}]);
            AreaNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.Area(perStageData{nfile,nch,nstage}.PeakAmp<0);AreaNeg{nfile,nch,nstage}]);
            AverageAmplitudeNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.AvgAmp(perStageData{nfile,nch,nstage}.PeakAmp<0);AverageAmplitudeNeg{nfile,nch,nstage}]);
            UpSlopeNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.UpSlope(perStageData{nfile,nch,nstage}.PeakAmp<0);UpSlopeNeg{nfile,nch,nstage}]);
            DownSlopeNeg{nfile,nch,nstage} = mean([perStageData{nfile,nch,nstage}.DownSlope(perStageData{nfile,nch,nstage}.PeakAmp<0);DownSlopeNeg{nfile,nch,nstage}]);
        end
    end
end

% Put it all in tables & write to excel
for nch = 1:length(PARAM.channels)
    for nstage = 1:length(PARAM.stages)
        % create tables
        stageNames = PARAM.stages{nstage};
        stageNames(isspace(stageNames)) = [];
        chNames = char(PARAM.channels(nch));
        chNames(strfind(chNames,'-')) = [];
        % all slow waves
        SummaryAll = table(ID', [NumberAllType{:,nch,nstage}]',[DurationAllType{:,nch,nstage}]',[FrequencyAllType{:,nch,nstage}]',[PeakAmplitudeAllType{:,nch,nstage}]',[AreaAllType{:,nch,nstage}]',[AverageAmplitudeAllType{:,nch,nstage}]',[UpSlopeAllType{:,nch,nstage}]',[DownSlopeAllType{:,nch,nstage}]');
        SummaryAll.Properties.VariableNames = {'ID',[chNames '_' stageNames '_Number'],[chNames '_' stageNames '_Duration'],[chNames '_' stageNames '_Frequency'],[chNames '_' stageNames '_PeakAmplitude'],[chNames '_' stageNames '_Area'],[chNames '_' stageNames '_AvgAmplitude'],[chNames '_' stageNames '_UpSlope'],[chNames '_' stageNames '_DownSlope']};
        % positive half waves
        SummaryPos = table(ID', [NumberPos{:,nch,nstage}]',[DurationPos{:,nch,nstage}]',[FrequencyPos{:,nch,nstage}]',[PeakAmplitudePos{:,nch,nstage}]',[AreaPos{:,nch,nstage}]',[AverageAmplitudePos{:,nch,nstage}]',[UpSlopePos{:,nch,nstage}]',[DownSlopePos{:,nch,nstage}]');
        SummaryPos.Properties.VariableNames = {'ID',[chNames '_' stageNames '_Number'],[chNames '_' stageNames '_Duration'],[chNames '_' stageNames '_Frequency'],[chNames '_' stageNames '_PeakAmplitude'],[chNames '_' stageNames '_Area'],[chNames '_' stageNames '_AvgAmplitude'],[chNames '_' stageNames '_UpSlope'],[chNames '_' stageNames '_DownSlope']};
        % negative half waves
        SummaryNeg = table(ID', [NumberNeg{:,nch,nstage}]',[DurationNeg{:,nch,nstage}]',[FrequencyNeg{:,nch,nstage}]',[PeakAmplitudeNeg{:,nch,nstage}]',[AreaNeg{:,nch,nstage}]',[AverageAmplitudeNeg{:,nch,nstage}]',[UpSlopeNeg{:,nch,nstage}]',[DownSlopeNeg{:,nch,nstage}]');
        SummaryNeg.Properties.VariableNames = {'ID',[chNames '_' stageNames '_Number'],[chNames '_' stageNames '_Duration'],[chNames '_' stageNames '_Frequency'],[chNames '_' stageNames '_PeakAmplitude'],[chNames '_' stageNames '_Area'],[chNames '_' stageNames '_AvgAmplitude'],[chNames '_' stageNames '_UpSlope'],[chNames '_' stageNames '_DownSlope']};
        % write to xlsx
        writetable(SummaryAll,[PARAM.resultDir filesep 'SlowWaveSummaryData_' chNames '_' char(PARAM.stages(nstage)) '.xlsx'],'Sheet','SummaryAll')
        writetable(SummaryPos,[PARAM.resultDir filesep 'SlowWaveSummaryData_' chNames '_' char(PARAM.stages(nstage)) '.xlsx'],'Sheet','SummaryPositive')
        writetable(SummaryNeg,[PARAM.resultDir filesep 'SlowWaveSummaryData_' chNames '_' char(PARAM.stages(nstage)) '.xlsx'],'Sheet','SummaryNegative')
    end
end
disp('ALL DONE!!!')

end
