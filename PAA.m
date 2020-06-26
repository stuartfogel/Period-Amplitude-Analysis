function [EEG] = PAA(EEG,ChOI,badData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function performs Period Amplitude Analysis (PAA) on EEGlab data to
% detect and measure slow waves. The events are marked on the EEG, saved
% as a new dataset and a results *.csv and *.mat tables are created.
%
% Input:
% ChOI = list of channel labels to run PAA, e.g., {'F3','Fz','F4'}
% badData = label for bad data event type, e.g., 'Movement'
%
% Requires:
% eeglab *.set file, ideally already sleep stage scored and movement
% artifacted (with events in the EEG.event structure).
%
% Output results table variables:
% 'ID': filename
% 'Channel': EEG channel (user-defined)
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
% June 24, 2020
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

%% INPUT ARGUMENTS
if nargin < 1; EEG = eeg_emptyset(); end
if nargin < 2; ChOI = 'F3 Fz F4'; end % default
if nargin < 3; badData = 'Movement'; end % name for movement artifact. Default: 'Movement'.

%% OPTIONAL PARAMETERS
eventName = {'SWpos','SWneg'}; % name of event. Default: {'SWpos','SWneg'}.
peaks = []; % mark event latencies at the HW peaks: 1, or at the zero-crossings: [] (default); note: marking peaks results in misalignment from results in output table. Useful for event-related analyses.

%% FILTER SETTING (default for slow wave detection 0.5-2Hz)
% DO NOT MODIFY THESE SETTINGS UNLESS YOU REALLY KNOW WHAT YOU ARE DOING!
% bandpass filter signal (filters with zero or linear phase-shifts needed to avoid signal distortions)
% e.g., frequency range of interest from 0.5 to 2 Hz
% band-pass filtered:
% sixth-order Chebyshev type II low-pass filter, -10 dB at 2.3 Hz (note: increased from 3 to 10 dB attenuation)
% third-order Chebyshev type II high-pass filter; -10 dB at 0.4 Hz (note: increased from 3 to 10 dB attenuation)
% The filters were applied in the forward and reverse directions, in order to achieve zero-phase distortion resulting in doubling of the filter order.
% The cut-off frequencies were selected to achieve minimal attenuation in the band of interest (0.5-2 Hz) and good attenuation at neighbouring frequencies, i.e. below 0.5 and above 2 Hz
% Based on doi: https://doi.org/10.1111/j.1365-2869.2009.00775.x
LPorder = 6;
LPfreq = 2.3;
LPattenuation = 10;
HPorder = 3;
HPfreq = 0.4;
HPattenuation = 10;
targetLPfreq = 2; % target lowpass frequency

%% PROCESS FILE
disp(['Processing dataset: ' EEG.setname '...'])

%% build the results table
ID = [];
CH = [];
firstHWlatency = [];
secondHWlatency = [];
segUpHWperiod = [];
segUpHWfreq = [];
segUpHWpeak = [];
segUpHWintamp = [];
segUpHWrecamp = [];
segUpHWupslope = [];
segUpHWdownslope = [];
segDownHWperiod = [];
segDownHWfreq = [];
segDownHWpeak = [];
segDownHWintamp = [];
segDownHWrecamp = [];
segDownHWupslope = [];
segDownHWdownslope = [];

SW = table(ID,CH,firstHWlatency,segUpHWperiod,segUpHWfreq,segUpHWpeak,segUpHWintamp,segUpHWrecamp,segUpHWupslope,segUpHWdownslope,secondHWlatency,segDownHWperiod,segDownHWfreq,segDownHWpeak,segDownHWintamp,segDownHWrecamp,segDownHWupslope,segDownHWdownslope,...
    'VariableNames',{'ID','Channel','firstLatency','posDur','posFreq','posPeakAmp','posArea','posAvgAmp','posUpSlope','posDownSlope','secondLatency','negDur','negFreq','negPeakAmp','negArea','negAvgAmp','negUpSlope','negDownSlope'});

%% find channels of interest
ChName = {EEG.chanlocs.labels};
nChOI = false(size(ChName));
ChOI = strsplit(ChOI,' ');
for n = 1:length(ChOI)
    nChOI = logical(nChOI + strcmp(ChName,ChOI{n}));
end
if ~any(nChOI)
    error('Missing or incorrect channel label input')
end
data = double(EEG.data(nChOI,:)); % filtfilt needs double precision data
clear ChName iChOI nChOI n

%% filter the EEG channels of interest
lp = designfilt('lowpassiir', 'FilterOrder', LPorder, 'StopbandFrequency', LPfreq, 'StopbandAttenuation', LPattenuation, 'SampleRate', EEG.srate);
hp = designfilt('highpassiir', 'FilterOrder', HPorder, 'StopbandFrequency', HPfreq, 'StopbandAttenuation', HPattenuation, 'SampleRate', EEG.srate);
datafilt = zeros(size(data));
for n=1:size(data,1)
    filtCh = filtfilt(lp,data(n,:));
    filtCh = filtfilt(hp,filtCh);
    datafilt(n,:) = filtCh(1,:);
end
clear hp lp LPorder LPfreq LPattenuation HPorder HPfreq HPattenuation filtCh n data

%% process each channel
for nch=1:size(datafilt,1)
    
    disp(['Processing channel: ' ChOI{nch} '...'])
    
    % find zero-crossings
    xdiff = diff(sign(datafilt(nch,:)));
    idx_up = find(xdiff>0);
    idx_down = find(xdiff<0);
    startDelay = min([idx_up(1),idx_down(1)]); % find the delay/onset of the first HW
    clear xdiff
    
    % segment data into half waves
    if length(idx_up) < length(idx_down) % so that you don't exceed the length of idx
        for n=1:length(idx_up)-1
            if idx_up(1) < idx_down(1) % which one is first?
                segUp{n} = datafilt(nch,idx_up(n)+1:idx_down(n));
                segDown{n} = datafilt(nch,idx_down(n)+1:idx_up(n+1));
            elseif idx_up(1) > idx_down(1)
                segUp{n} = datafilt(nch,idx_down(n)+1:idx_up(n));
                segDown{n} = datafilt(nch,idx_up(n)+1:idx_down(n+1));
            end
        end
    elseif length(idx_up) > length(idx_down) % so that you don't exceed the length of idx
        for n=1:length(idx_down)-1
            if idx_up(1) < idx_down(1) % which one is first?
                segUp{n} = datafilt(nch,idx_up(n)+1:idx_down(n));
                segDown{n} = datafilt(nch,idx_down(n)+1:idx_up(n+1));
            elseif idx_up(1) > idx_down(1)
                segUp{n} = datafilt(nch,idx_down(n)+1:idx_up(n));
                segDown{n} = datafilt(nch,idx_up(n)+1:idx_down(n+1));
            end
        end
    end
    
    % identify and count HWs > peak amplitude threshold & longer than high freq cut-off
    % e.g., slow waves: 75 uV peak-to-peak, 37.5 uV for positive or
    % negative half-waves e.g., > 0.25 sec (re; impossible to have a SW <2Hz happen in that time)
    segUpDetect = false(1,length(segUp{n}));
    for n=1:length(segUp)
        segUpDetect(n) = logical(max(segUp{n})>37.5 && length(segUp{n})>1/targetLPfreq*EEG.srate/2 && min(segDown{n})<-37.5 && length(segDown{n})>1/targetLPfreq*EEG.srate/2);
    end
    segDownDetect = false(1,length(segDown{n}));
    for n=1:length(segDown)
        segDownDetect(n) = logical(min(segDown{n})<-37.5 && length(segDown{n})>1/targetLPfreq*EEG.srate/2 && max(segUp{n})>37.5 && length(segUp{n})>1/targetLPfreq*EEG.srate/2);
    end
    
    % compute HW period, frequency, peak amplitude, peak, integrated (i.e., area under curve: sum of rectified values * 1/srate),
    % rectified amplitude (i.e., avg amplitude: integrated amplitude/HW period) and up (crossing to peak) and down (peak to crossing) slope
    
    segUpHWperiod = NaN(size(segUp));
    segUpHWfreq = NaN(size(segUp));
    segUpHWpeak = NaN(size(segUp));
    segUpHWintamp = NaN(size(segUp));
    segUpHWrecamp = NaN(size(segUp));
    segUpHWpeakLat = NaN(size(segUp));
    segUpHWupslope = NaN(size(segUp));
    segUpHWdownslope = NaN(size(segUp));
    
    for n=1:length(segUp)
        if segUpDetect(n) == 1
            segUpHWperiod(n) = length(segUp{n})/EEG.srate; % compute HW period
            segUpHWfreq(n) = 1/segUpHWperiod(n); % compute frequency of HW (1/period)
            segUpHWpeak(n) = max(segUp{n}); % compute peak amplitude
            segUpHWintamp(n) = abs(sum(segUp{n}))*1/EEG.srate; % compute integrated amplitude
            segUpHWrecamp(n) = segUpHWintamp(n)/segUpHWperiod(n); % compute rectified amplitude
            segUpHWpeakLat(n) = find(segUp{n} == max(segUp{n})); % compute peak latency
            segUpHWupslope(n) = segUpHWpeak(n)/(segUpHWpeakLat(n)/EEG.srate*1000); % compute upward slope in uV/ms
            segUpHWdownslope(n) = segUpHWpeak(n)/(length(segUp{n})-segUpHWpeakLat(n)/EEG.srate*1000); % compute downward slope in uV/ms
        end
    end
    
    segDownHWperiod = NaN(size(segDown));
    segDownHWfreq = NaN(size(segDown));
    segDownHWpeak = NaN(size(segDown));
    segDownHWintamp = NaN(size(segDown));
    segDownHWrecamp = NaN(size(segDown));
    segDownHWpeakLat = NaN(size(segUp));
    segDownHWupslope = NaN(size(segDown));
    segDownHWdownslope = NaN(size(segDown));
    
    for n=1:length(segDown)
        if segDownDetect(n) == 1
            segDownHWperiod(n) = length(segDown{n})/EEG.srate; % compute HW period
            segDownHWfreq(n) = 1/segDownHWperiod(n); % compute frequency of HW (1/period)
            segDownHWpeak(n) = min(segDown{n}); % compute peak amplitude
            segDownHWintamp(n) = abs(sum(segDown{n}))*1/EEG.srate; % compute integrated amplitude
            segDownHWrecamp(n) = segDownHWintamp(n)/segDownHWperiod(n); % compute rectified amplitude
            segDownHWpeakLat(n) = find(segDown{n} == min(segDown{n})); % compute peak latency
            segDownHWupslope(n) = segDownHWpeak(n)/(segDownHWpeakLat(n)/EEG.srate*1000); % compute upward slope in uV/ms
            segDownHWdownslope(n) = segDownHWpeak(n)/(length(segDown{n})-segDownHWpeakLat(n)/EEG.srate*1000); % compute downward slope in uV/ms
        end
    end
    
    %% mark events on the EEG
    if ~isfield(EEG.event,'channel')
        for i=1:length(EEG.event)
            EEG.event(i).channel = ''; % add empty channel field
        end
    end
    % create empty EEGlab event structure to match existing to merge
    fields = fieldnames(EEG.event)';
    fields{2,1} = {};
    eventsUp = struct(fields{:});
    eventsDown = struct(fields{:});
    clear fields
    
    % add each HW SW event to the empty EEGlab event structure
    if idx_up(1) < idx_down(1) % which one is first?
        latency = startDelay;
        nevt = 1; % valid event counter
        for n = 1:length(segUpDetect)
            if segUpDetect(n) == 1
                eventsUp(nevt).type = eventName{1};
                if isempty(peaks)
                    eventsUp(nevt).latency = latency;
                    eventsUp(nevt).duration = segUpHWperiod(n)*EEG.srate;
                else
                    eventsUp(nevt).latency = latency + segUpHWpeakLat(n);
                    eventsUp(nevt).duration = 1; % one sample peak marker
                end
                eventsUp(nevt).channel = ChOI{nch};
                eventsUp(nevt).urevent = [];
                nevt = nevt + 1; % valid event counter
            end
            latency = latency + length(segUp{n}) + length(segDown{n});
        end
        latency = startDelay + length(segUp{1});
        nevt = 1; % valid event counter
        for n = 1:length(segDownDetect)
            if segDownDetect(n) == 1
                eventsDown(nevt).type = eventName{2};
                if isempty(peaks)
                    eventsDown(nevt).latency = latency;
                    eventsDown(nevt).duration = segDownHWperiod(n)*EEG.srate;
                else
                    eventsDown(nevt).latency = latency + segDownHWpeakLat(n);
                    eventsDown(nevt).duration = 1; % one sample peak marker
                end
                eventsDown(nevt).channel = ChOI{nch};
                eventsDown(nevt).urevent = [];
                nevt = nevt + 1; % valid event counter
            end
            if n == length(segUp) % catch exceed end of array
                break
            else
                latency = latency + length(segUp{n+1}) + length(segDown{n});
            end
        end
        clear latency
    elseif idx_up(1) > idx_down(1)
        latency = startDelay + length(segDown{1});
        nevt = 1; % valid event counter
        for n = 1:length(segUpDetect)
            if segUpDetect(n) == 1
                eventsUp(nevt).type = eventName{1};
                if isempty(peaks)
                    eventsUp(nevt).latency = latency;
                    eventsUp(nevt).duration = segUpHWperiod(n)*EEG.srate;
                else
                    eventsUp(nevt).latency = latency + segUpHWpeakLat(n);
                    eventsUp(nevt).duration = 1; % one sample peak marker
                end
                eventsUp(nevt).channel = ChOI{nch};
                eventsUp(nevt).urevent = [];
                nevt = nevt + 1; % valid event counter
            end
            if n == length(segDown) % catch exceed end of array
                break
            else
                latency = latency + length(segUp{n}) + length(segDown{n+1});
            end
        end
        latency = startDelay;
        nevt = 1; % valid event counter
        for n = 1:length(segDownDetect)
            if segDownDetect(n) == 1
                eventsDown(nevt).type = eventName{2};
                if isempty(peaks)
                    eventsDown(nevt).latency = latency;
                    eventsDown(nevt).duration = segDownHWperiod(n)*EEG.srate;
                else
                    eventsDown(nevt).latency = latency + segDownHWpeakLat(n);
                    eventsDown(nevt).duration = 1; % one sample peak marker
                end
                eventsDown(nevt).channel = ChOI{nch};
                eventsDown(nevt).urevent = [];
                nevt = nevt + 1; % valid event counter
            end
            latency = latency + length(segUp{n}) + length(segDown{n});
        end
        clear latency
    end
    
    % concatenate event structures and sort
    EEG.event = [EEG.event eventsUp eventsDown];
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
    
    %% populate table with results
    ID = [];
    CH = [];
    nevt = 1; % valid event counter
    if length(segDownDetect) > length(segUpDetect) % just in case they differ in length, which should never be the case!
        for n=1:length(segDownDetect)
            ID{n} = EEG.setname;
            CH{n} = ChOI{nch};
            if idx_up(1) < idx_down(1) % which one is first?
                if segUpDetect(n) == 1
                    firstHWlatency{n} = eventsUp(nevt).latency;
                    secondHWlatency{n} = eventsDown(nevt).latency;
                    nevt = nevt + 1; % valid event counter
                else
                    firstHWlatency{n} = NaN;
                end
            elseif idx_up(1) > idx_down(1)
                if segDownDetect(n) == 1
                    firstHWlatency{n} = eventsDown(nevt).latency;
                    secondHWlatency{n} = eventsUp(nevt).latency;
                    nevt = nevt + 1; % valid event counter
                else
                    firstHWlatency{n} = NaN;
                    secondHWlatency{n} = NaN;
                end
            end
        end
    else
        for n=1:length(segUpDetect)
            ID{n} = EEG.setname;
            CH{n} = ChOI{nch};
            if idx_up(1) < idx_down(1) % which one is first?
                if segUpDetect(n) == 1
                    firstHWlatency{n} = eventsUp(nevt).latency;
                    secondHWlatency{n} = eventsDown(nevt).latency;
                    nevt = nevt + 1; % valid event counter
                else
                    firstHWlatency{n} = NaN;
                    secondHWlatency{n} = NaN;
                end
            elseif idx_up(1) > idx_down(1)
                if segDownDetect(n) == 1
                    firstHWlatency{n} = eventsDown(nevt).latency;
                    secondHWlatency{n} = eventsUp(nevt).latency;
                    nevt = nevt + 1; % valid event counter
                else
                    firstHWlatency{n} = NaN;
                    secondHWlatency{n} = NaN;
                end
            end
        end
    end
    SWnew = table(ID',CH',firstHWlatency',segUpHWperiod',segUpHWfreq',segUpHWpeak',segUpHWintamp',segUpHWrecamp',segUpHWupslope',segUpHWdownslope',secondHWlatency',segDownHWperiod',segDownHWfreq',segDownHWpeak',segDownHWintamp',segDownHWrecamp',segDownHWupslope',segDownHWdownslope',...
        'VariableNames',{'ID','Channel','firstLatency','posDur','posFreq','posPeakAmp','posArea','posAvgAmp','posUpSlope','posDownSlope','secondLatency','negDur','negFreq','negPeakAmp','negArea','negAvgAmp','negUpSlope','negDownSlope'});
    SW = [SW ; SWnew];
    
    clear segUpDetect segDownDetect segUp segDown idx_up idx_down firstHWlatency secondHWlatency
    clear segUpHWperiod segUpHWfreq segUpHWpeak segUpHWintamp segUpHWrecamp segUpHWupslope segUpHWdownslope segDownHWperiod segDownHWfreq segDownHWpeak segDownHWintamp segDownHWrecamp segDownHWupslope segDownHWdownslope segDownHWpeakLat segUpHWpeakLat
    clear eventsUp eventsDown i n nevt
    
end

%% mark for delete SW during movement artifact
Event = EEG.event;
evtIdx = find(ismember({Event.type},eventName));
ToRmv = [];

for iEvt = evtIdx % loop on event
    bad = find(ismember({Event(1:iEvt).type},badData),1,'last');
    if ~isempty(bad)
        if Event(bad).latency + Event(bad).duration > Event(iEvt).latency
            ToRmv(end+1) = iEvt;
        end
    end
end

% let's flag the bad event
for iEvt = ToRmv
    Event(iEvt).type = 'BAD'; % so that they will get flagged for removal from results table and EEG.event below
end
EEG.event = Event;
clear bad Event

% remove NaN values from the results table
SW=SW(~any(ismissing(SW),2),:);

% find the SW marked as bad in the table by their latency
iRow=[];
for nEvt=ToRmv
    row = find([SW.firstLatency{:}]' == EEG.event(nEvt).latency);
    if ~isempty(row)
        iRow(end+1) = row;
    end
end
% remove the bad events from the results table
SW(iRow,:) = [];

% now, lets remove the bad events from EEG.event
EEG.event(ToRmv) = []; % delete

% check that everything is ok
EEG = eeg_checkset(EEG, 'eventconsistency');
eeg_checkset(EEG);


%% SELECT OUTPUT DIRECTORY
disp('Please select a directory in which to save the results.');
resultDir = uigetdir('', 'Select the directory in which to save the results');

%% save out the results
writetable(SW,[resultDir filesep 'SWevents.csv'],'Delimiter',',');
save([resultDir filesep 'SWevents.mat'], 'SW');
clear SW SWnew startDelay nch datafilt ToRmv Event

end