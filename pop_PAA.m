function [EEG,com] = pop_PAA(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function performs Period Amplitude Analysis (PAA) on EEGlab data to
% detect and measure slow waves. The events are marked on the EEG, saved
% as a new dataset and a results *.csv and *.mat tables are created.
%
% Input:
% ChOI = list of channel labels to run PAA, e.g., {'F3','Fz','F4'}
% badData = label for bad data event type, e.g., 'Movement'
% allSleepStages = all sleep stages included in scoring, e.g., {'N1','N2','N3','R','W'}
%
% Requires:
% eeglab *.set file, ideally already sleep stage scored and movement
% artifacted (with events in the EEG.event structure).
%
% Output results table variables:
% 'N': file number
% 'ID': filename
% 'Channel': EEG channel (user-defined)
% 'sleepStage': Sleep Stage that the first part of the half wave occurs in
% 'Latency': latency of first half-wave in data points
% 'Duration': positive half-wave period in sec
% 'Frequency': positive half-wave ocsillatory frequency in Hz
% 'PeakAmp': positive half-wave peak amplitude in uV
% 'Area': positive half-wave area under the curve (i.e., integrated amplitude = sum of rectified values * 1/srate) in uV*sec
% 'AvgAmp': positive half-wave average amplitude (i.e., rectified amplitude = integrated amplitude/period) in uV
% 'UpSlope': positive half-wave upward slope in uV/sec
% 'DownSlope': positive half-wave downward slope in uV/sec
%
% June 24, 2020 Version 1.0
% Aug  27, 2020 Revised 1.1 Critical bug fixes: ch order, polarity, channel
%   labels
% Sept 13, 2020 Revised 1.2 included sleep stages in output and SW events,
%   fixed bug for SW inclusion criteria, optimised code
% Sept 16, 2020 Revised 1.3 negative slope calculation bug fixed. Improved
%   detection criteria to include any adjacent HWs
% Sept 23, 2020 Revised 1.4 major fix for starting issue with polarity and
%   table creation for multiple files.
% Sept 24, 2020 Revised 1.5 fixed conflict with identical latency events
% Sept 28, 2020 Revised 1.6 adjusted filtering parameters and functions to
%   improve filter response - AG
% Jun 1, 2021 Revised 1.7 (incorporated updates from AG): 
%   1. updates to HW threshold process. The peak-to-peak amplitude is now 
%   calculated by finding the difference between the maximum absolute peak 
%   of each HW and it's oppositely valanced HW neighbours. As before, the 
%   length of each HW and it's neighbour is also checked to ensure it falls
%   within the correct frequency range. Only HWs that meet both criteria 
%   are included in future analyses. 
%   2. added feature to save each subjects data as separate .csv files in 
%   addition to concatinated .csv for all subjects 
%   3. minor fix to correct mismatch between filename and setname in output 
%   table - SF
%   4. added half wave amplitude threshold in addition to p2p amplitude 
%   threshold. - AG
%   5. added functionality to remove unwanted SW during sleepstages of 
%   non-interest - SF
% Jun 21, 2021 Revised 1.8: added feature to remove SW events outside
%   Lights OFF/ON tags.
%
% Copyright, Sleep Well. https://www.sleepwellpsg.com
%
% Method adpated from:
%
% Feinberg et al, 1978. https://doi.org/10.1016/0013-4694(78)90266-3
% Geering et al, 1993. https://doi.org/10.1111/j.1365-2869.1993.tb00074.x
% Bersagliere and Achermann, 2010. https://doi.org/10.1111/j.1365-2869.2009.00775.x
% ...also, see recent papers by Tononi's group
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle history
com = '';
if nargin < 1
    help pop_PAA;
    return;
end

% GUI geometry setup
g = [3, 2];
geometry = {1,g,g,g,g,1,[2 2 1]};
geomvert = [1 1 1 1 1 1 1];

% select channels
cb_chan = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''ChOI''));';

% build GIU
uilist = { ...
    ... label settings
    {'style', 'text', 'string', 'PAA for slow wave (0.5-2Hz, >75uV) detection'} ...
    {'style', 'text', 'string', 'Label for all sleep stages'} ...
    {'style', 'edit', 'string', 'N1 N2 N3 REM Wake' 'tag' 'allSleepStages'} ...
    {'style', 'text', 'string', 'Label for sleep stages to exclude SW'} ...
    {'style', 'edit', 'string', 'N1 REM Wake' 'tag' 'badSleepstages'} ...
    {'style', 'text', 'string', 'Label for bad data'} ...
    {'style', 'edit', 'string', 'Movement' 'tag' 'badData'} ...
    {'style', 'text', 'string', 'Labels for lights ON/OFF tags (comma separated)'} ...
    {'style', 'edit', 'string', 'Lights Off, Lights On' 'tag' 'lightsTags'} ...
    ... channel options
    { 'style' 'text'       'string' '' } ...
    { 'style' 'text'       'string' 'Channel labels or indices' } ...
    { 'style' 'edit'       'string' 'F3 Fz F4' 'tag' 'ChOI' }  ...
    { 'style' 'pushbutton' 'string' '...' 'callback' cb_chan }
    };

% channel labels
% --------------
if ~isempty(EEG(1).chanlocs)
    tmpchanlocs = EEG(1).chanlocs;
else
    tmpchanlocs = [];
    for index = 1:EEG(1).nbchan
        tmpchanlocs(index).labels = int2str(index);
        tmpchanlocs(index).type = '';
    end
end

% launch gui
result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Period Amplitude Analysis -- PAA', 'helpcom', 'pophelp(''PAA'')', 'userdata', tmpchanlocs);

% launch PAA
if ~isempty(result)
    allSleepStages = result{1};
    badSleepstages = result{2};    
    badData = result{3};
    lightsTags = result{4};
    ChOI = result{5};
    % launch pipeline
    [EEG] = PAA(EEG,ChOI,badData,allSleepStages,badSleepstages,lightsTags);
else
    com = '';
    return
end

com = sprintf('EEG = PAA(%s,%s,%s,%s,%s,%s);',inputname(1),ChOI,badData,allSleepStages,badSleepstages,lightsTags);

end