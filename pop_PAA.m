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
% 'ID': filename
% 'Channel': EEG channel (user-defined)
% 'sleepStage': Sleep Stage that the first part of the half wave occurs in
% 'firstLatency': latency of first half-wave in data points
% 'posDur': positive half-wave period in sec
% 'posFreq': positive half-wave ocsillatory frequency in Hz
% 'posPeakAmp': positive half-wave peak amplitude in uV
% 'posArea': positive half-wave area under the curve (i.e., integrated amplitude = sum of rectified values * 1/srate) in uV*sec
% 'posAvgAmp': positive half-wave average amplitude (i.e., rectified amplitude = integrated amplitude/period) in uV
% 'posUpSlope': positive half-wave upward slope in uV/sec
% 'posDownSlope': positive half-wave downward slope in uV/sec
% 'secondLatency': latency of second half-wave in data points
% 'negDur': negative half-wave period in sec
% 'negFreq': negative half-wave ocsillatory frequency in Hz
% 'negPeakAmp': negative half-wave peak amplitude in uV
% 'negArea': negative half-wave area under the curve (i.e., integrated amplitude = sum of rectified values * 1/srate) in uV*sec
% 'negAvgAmp': negative half-wave average amplitude (i.e., rectified amplitude = integrated amplitude/period) in uV
% 'negUpSlope': negative half-wave upward slope in uV/sec
% 'negDownSlope': negative half-wave downward slope in uV/sec
%
% June 24, 2020 Version 1.0
% Aug 27, 2020  Revised 1.1 Critical bug fixes: ch order, polarity, channel
% labels
% Sept 13, 2020  Revised 1.2 included sleep stages in output and SW events,
% fixed bug for SW inclusion criteria, optimised code
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
geometry = {1,g,g,1,[2 2 1]};
geomvert = [1 1 1 1 1];

% select channels
cb_chan = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''ChOI''));';

% build GIU
uilist = { ...
    ... label settings
    {'style', 'text', 'string', 'PAA for slow wave (0.5-2Hz, >75uV) detection'} ...
    {'style', 'text', 'string', 'Label for sleep stages'} ...
    {'style', 'edit', 'string', 'N1 N2 N3 REM Wake' 'tag' 'allSleepStages'} ...
    {'style', 'text', 'string', 'Label for bad data'} ...
    {'style', 'edit', 'string', 'Movement' 'tag' 'badData'} ...
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
    badData = result{2};
    ChOI = result{3};
    % launch pipeline
    [EEG] = PAA(EEG,ChOI,badData,allSleepStages);
else
    com = '';
    return
end

com = sprintf('EEG = PAA(%s,%s,%s,%s);',inputname(1),ChOI,badData,allSleepStages);

end