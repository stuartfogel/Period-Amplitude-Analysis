function [EEG] = PAA_standalone(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function performs Period Amplitude Analysis (PAA) on EEGlab data to
% detect and measure slow waves. The events are marked on the EEG, saved
% as a new dataset and a results *.csv and *.mat tables are created.
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
% Aug  27, 2020 Revised 1.1: Critical bug fixes: ch order, polarity, channel
%   labels
% Sept 13, 2020 Revised 1.2: included sleep stages in output and SW events,
%   fixed bug for SW inclusion criteria, optimised code
% Sept 16, 2020 Revised 1.3: negative slope calculation bug fixed. Improved
%   detection criteria to include any adjacent HWs
% Sept 23, 2020 Revised 1.4: major fix for starting issue with polarity and
%   table creation for multiple files.
% Sept 24, 2020 Revised 1.5: fixed conflict with identical latency events
% Sept 28, 2020 Revised 1.6: adjusted filtering parameters and functions to
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

%% LOAD EEGLAB
eeglab % make sure all eeglab functions are added to the path
clear; close all; clc;

%% USER-DEFINED PARAMETERS
ChOI = {'Fz'}; % user-defined. Default: {'Fz','Cz','Pz'}.
peaks = []; % mark event latencies at the HW peaks 1, or at the zero-crossings [] (default); note: marking peaks results in misalignment from results in output table
eventName = {'SWpos','SWneg'}; % name of event. Default: {'SWpos','SWneg'}.
allSleepStages = {'N1','N2','N3','R','W'}; % all sleep stages included in scoring. Default: {'N1','N2','N3','REM','Wake'}.
badSleepstages = {'N1','R','W'}; % sleep stages excluded from SW detection, e.g., {'N1','REM','Wake'}.
badData = 'Movement'; % name for movement artifact. Default: 'Movement'.
lightsTags = {'Lights Off','Lights On'}; % tags for lights on and lights off. used to remove false detections outside lights off interval. Default: {'Lights Off','Lights On'}.

%% FILTER SETTING (default for slow wave detection 0.5-4Hz)
% DO NOT MODIFY THESE SETTINGS UNLESS YOU REALLY KNOW WHAT YOU ARE DOING!
% bandpass filter signal (filters with zero or linear phase-shifts needed to avoid signal distortions)
% e.g., frequency range of interest from 0.5 to 4 Hz
% band-pass filtered:
% 32nd-order Chebyshev type II low-pass filter, -80 dB stopband attenuation (note: increased from 10 to 80 dB attenuation)
% 64th-order Chebyshev type II high-pass filter; -80 dB stopband attenuation (note: increased from 10 to 80 dB attenuation)
% The filters were applied in the forward and reverse directions, to achieve zero-phase distortion resulting in doubling of the filter order.
% The cut-off frequencies were selected to achieve minimal attenuation in the band of interest (0.5-4 Hz) and good attenuation at neighbouring frequencies, i.e. below 0.5 and above 4 Hz
% Adapted from doi: https://doi.org/10.1111/j.1365-2869.2009.00775.x
LPorder  = 32;   % lowpass filter order
LPfreq   = 2.15; % use 2.15 for target of 2 Hz, or 4.3 for target of 4Hz
HPorder  = 64;   % highpass filter order
HPfreq   = 0.46; % highpass filter cutoff
filtAttn = 80;   % stopband attenuation (dB)

targetLPfreq = 2; % target lowpass frequency: use 2 for target of 2 Hz, or use 4 for target of 4Hz
% use fvtool(lp,hp) after building filters below to visualize

%% MANUALLY SELECT *.SET/MAT FILES
disp('Please select file(s) to process.');
[filename, pathname] = uigetfile2( ...
    {'*.set','EEGlab Dataset (*.set)'; ...
     '*.mat','MAT-files (*.mat)'; ...
     '*.*',  'All Files (*.*)'}, ...
     'multiselect', 'on');

% check the filename(s)
if isequal(filename,0) || isequal(pathname,0) % no files were selected
    disp('User selected Cancel')
    return;
else
    if ischar(filename) % Only one file was selected.
        filename = cellstr(filename); % Put the filename in the same cell structure as multiselect
    end
end

%% SELECT OUTPUT DIRECTORY
disp('Please select a directory in which to save the results.');
resultDir = uigetdir('', 'Select the directory in which to save the results');

tic

%% build the results table
N = [];
ID = [];
CH = [];
sleepStage = [];
HWlatency = [];
HWperiod = [];
HWfreq = [];
HWpeak = [];
HWintamp = [];
HWrecamp = [];
HWupslope = [];
HWdownslope = [];

SW = table(N,ID,CH,sleepStage,HWlatency,HWperiod,HWfreq,HWpeak,HWintamp,HWrecamp,HWupslope,HWdownslope,...
    'VariableNames',{'N','ID','Channel','sleepStage','Latency','Duration','Frequency','PeakAmp','Area','AvgAmp','UpSlope','DownSlope'});

%% PROCESS ALL FILES
for nfile = 1:length(filename)
    
    disp(['Processing file: ' filename{nfile} '...'])
    
    %% load the EEGlab dataset
    EEG = pop_loadset('filename',filename{1,nfile},'filepath',pathname);
    
    %% find channels of interest
    ChName = {EEG.chanlocs.labels};
    nChOI = false(size(ChName));
    for n = 1:length(ChOI)
        chIdx(n) = find(strcmp(ChName, ChOI{n}));
        nChOI = logical(nChOI + strcmp(ChName,ChOI{n}));
    end
    if ~any(nChOI)
        error('Missing or incorrect channel label input')
    end
    data = double(EEG.data(chIdx,:)); % filtfilt needs double precision data. Note: this also re-orders channels to match ChOI channel order.
    clear n nChOI chIdx
    
    %% filter the EEG channels of interest
    disp('Filtering the data...')
    
    % build filters
    lp = design(fdesign.lowpass('N,Fst,Ast', LPorder, LPfreq, filtAttn, EEG.srate), 'cheby2');
    hp = design(fdesign.highpass('N,Fst,Ast', HPorder, HPfreq, filtAttn, EEG.srate), 'cheby2');
    % fvtool(lp,hp) % use to visualize filters
    
    % filter the data
    datafilt = zeros(size(data)); % create array for filtered data
    for n=1:size(data,1)
        filtCh = filtfilt(lp.sosMatrix,lp.ScaleValues, data(n,:)); % low-pass filter
        filtCh = filtfilt(hp.sosMatrix, hp.ScaleValues, filtCh); % high-pass filter
        datafilt(n,:) = filtCh(1,:); % add filtered data to array
    end
    clear hp lp filtCh n data
    
    %% process each channel
    for nch=1:size(datafilt,1)
        
        disp(['Processing channel: ' num2str(nch) ' of ' num2str(size(datafilt,1)) ' ...'])
        
        % find zero-crossings
        disp('Finding zero crossings...')
        xdiff = diff(sign(datafilt(nch,:)));
        idx_up = find(xdiff>0);
        idx_down = find(xdiff<0);
        startDelay = min([idx_up(1),idx_down(1)]); % find the delay/onset of the first half wave (HW)
        clear xdiff
        
        % segment data into half waves
        disp('Segmenting into half-waves...')
        for n=1:min(length(idx_up),length(idx_down))-1 % take min so that you don't exceed the length of idx
            if idx_up(1) < idx_down(1) % which one is first?
                segUp{n} = datafilt(nch,idx_up(n)+1:idx_down(n));
                segDown{n} = datafilt(nch,idx_down(n)+1:idx_up(n+1));
            elseif idx_up(1) > idx_down(1)
                segDown{n} = datafilt(nch,idx_down(n)+1:idx_up(n));
                segUp{n} = datafilt(nch,idx_up(n)+1:idx_down(n+1));
            end
        end
        clear n
        
        % verify that polarity of HWs are correct
        if sign(segUp{1}(1)) ~= 1 && sign(segDown{1}(1)) ~= -1
            for nSeg=1:length(segUp)
                segUp{nSeg}=segUp{nSeg}*-1;
                segDown{nSeg}=segDown{nSeg}*-1;
            end
        end
        
        % identify and count HWs > peak amplitude threshold & longer than high freq cut-off
        % e.g., slow waves: 75 uV peak-to-peak, 37.5 uV for positive or negative half-waves 
        % e.g., > 0.125 sec (re; impossible to have a SW <4Hz happen in that time)
        p2pMin = 75; 
        HWmin  = 37.5; 
        minLen = EEG.srate/(targetLPfreq.*2);  % minimum length of HW in samples
        
        % % % detect up segments
        upPeak = cellfun(@(p) max(abs(p)), segUp); % find abs max amplitude of up segments
        upLen  = cellfun(@length, segUp);   % length of up segments in samples
        
        dnPeak = cellfun(@(n) max(abs(n)), segDown); % find abs max amplitude of down segments
        dnLen  = cellfun(@length, segDown); % length of down segments in samples
        
        if idx_up(1) < idx_down(1)  % if first HW is postitive
            dnPeak = [nan, dnPeak(1:end-1); dnPeak];    % make array of peaks for negative HW before and after postitive wave
            dnLen  = [nan, dnLen(1:end-1); dnLen];      % make array of segment lengths for neg HW before and after positive HW
        else                        % if first HW is negative
            dnPeak = [dnPeak; dnPeak(2:end), nan];    % make array of peaks for negative HW before and after postitive wave
            dnLen  = [dnLen; dnLen(2:end), nan];      % make array of segment lengths for neg HW before and after positive HW
        end
        
        % identify positive HWs by amplitude
        HWamp  = sum([upPeak; max(dnPeak, [], 1)], 1) >= p2pMin & ... % add max abs amplitude of up segment to max abs amplitude of neighbouring down segs to determine p2p range
            all([upPeak; max(dnPeak, [], 1)] >= HWmin, 1);            % check that each HW segment is exceeds HW amplitude threshold
        
        % identify positive HWs by length
        HW_idx = dnPeak == repmat(max(dnPeak, [], 1), size(dnPeak, 1), 1); % index of down segments used in peak-to-peak calculation
        HWlen  = upLen >= minLen & dnLen(HW_idx)' >= minLen; 
        
        segUpDetect = HWamp & HWlen; % index of positive HWs that meet both amplitude and length thresholds
        clear upPeak upLen dnPeak dnLen HWamp HW_idx HWlen % clear variables so they don't interfere with next section
        
        % % % detect down segments
        dnPeak = cellfun(@(n) max(abs(n)), segDown);
        dnLen  = cellfun(@length, segDown); 
        
        upPeak = cellfun(@(p) max(abs(p)), segUp); 
        upLen  = cellfun(@length, segUp); 
        
        if idx_down(1) < idx_up(1)  % if first HW is negative
            upPeak = [nan, upPeak(1:end-1); upPeak];    % make array of peaks for pos HW before and after neg HW
            upLen  = [nan, upLen(1:end-1); upLen];      % make array of segment lengths for pos HW before and after neg HW
        else                        % if first HW is positive
            upPeak = [upPeak; upPeak(2:end), nan];    % make array of peaks for pos HW before and after neg wave
            upLen  = [upLen; upLen(2:end), nan];      % make array of segment lengths for pos HW before and after neg HW
        end
        
        % identify negative HWs by amplitude
        HWamp  = sum([dnPeak; max(upPeak, [], 1)], 1) >= p2pMin & ... % add max abs amplitude of neg HW to max abs amplitude of neighbouring up segs to determine p2p range
            all([dnPeak; max(upPeak, [], 1)] >= HWmin, 1);           % check that each HW segment is exceeds HW amplitude threshold
        
        % identify negative HWs by length
        HW_idx = upPeak == repmat(max(upPeak, [], 1), size(upPeak, 1), 1); 
        HWlen  = dnLen >= minLen & upLen(HW_idx)' >= minLen; 
        
        % save index of negative HWs that meet amplitude and length thresholds
        segDownDetect = HWamp & HWlen; 
        clear upPeak upLen dnPeak dnLen HWamp HW_idx HWlen % clear variables
                
        % compute HW period, frequency, peak amplitude, peak, integrated (i.e., area under curve: sum of rectified values * 1/srate),
        % rectified amplitude (i.e., avg amplitude: integrated amplitude/HW period) and up (crossing to peak) and down (peak to crossing) slope
        
        disp('Measuring half-waves...')
        
        upHWperiod = NaN(1,length(segUp));
        upHWfreq = NaN(1,length(segUp));
        upHWpeak = NaN(1,length(segUp));
        upHWintamp = NaN(1,length(segUp));
        upHWrecamp = NaN(1,length(segUp));
        upHWpeakLat = NaN(1,length(segUp));
        upHWupslope = NaN(1,length(segUp));
        upHWdownslope = NaN(1,length(segUp));
        
        for n=1:length(segUp)
            if segUpDetect(n) == 1
                upHWperiod(n) = length(segUp{n})/EEG.srate; % compute HW period
                upHWfreq(n) = 1/upHWperiod(n); % compute frequency of HW (1/period)
                upHWpeak(n) = max(segUp{n}); % compute peak amplitude
                upHWintamp(n) = abs(sum(segUp{n}))*1/EEG.srate; % compute integrated amplitude
                upHWrecamp(n) = upHWintamp(n)/upHWperiod(n); % compute rectified amplitude
                upHWpeakLat(n) = find(segUp{n} == max(segUp{n})); % compute peak latency
                upHWupslope(n) = upHWpeak(n)/(upHWpeakLat(n)/EEG.srate*1000); % compute upward slope in uV/ms
                upHWdownslope(n) = (upHWpeak(n)/((length(segUp{n})-upHWpeakLat(n))/EEG.srate*1000))*-1; % compute downward slope in uV/ms
            end
        end
        
        downHWperiod = NaN(1,length(segDown));
        downHWfreq = NaN(1,length(segDown));
        downHWpeak = NaN(1,length(segDown));
        downHWintamp = NaN(1,length(segDown));
        downHWrecamp = NaN(1,length(segDown));
        downHWpeakLat = NaN(1,length(segDown));
        downHWupslope = NaN(1,length(segDown));
        downHWdownslope = NaN(1,length(segDown));
        
        for n = 1:length(segDown)
            if segDownDetect(n) == 1
                downHWperiod(n) = length(segDown{n})/EEG.srate; % compute HW period
                downHWfreq(n) = 1/downHWperiod(n); % compute frequency of HW (1/period)
                downHWpeak(n) = min(segDown{n}); % compute peak amplitude
                downHWintamp(n) = abs(sum(segDown{n}))*1/EEG.srate; % compute integrated amplitude
                downHWrecamp(n) = downHWintamp(n)/downHWperiod(n); % compute rectified amplitude
                downHWpeakLat(n) = find(segDown{n} == min(segDown{n})); % compute peak latency
                downHWupslope(n) = (downHWpeak(n)/(downHWpeakLat(n)/EEG.srate*1000))*-1; % compute upward slope in uV/ms
                downHWdownslope(n) = downHWpeak(n)/((length(segDown{n})-downHWpeakLat(n))/EEG.srate*1000); % compute downward slope in uV/ms
            end
        end
        
        % concatenate up and down segments
        HWperiod = [upHWperiod, downHWperiod];
        HWfreq = [upHWfreq, downHWfreq];
        HWpeak = [upHWpeak, downHWpeak];
        HWintamp = [upHWintamp, downHWintamp];
        HWrecamp = [upHWrecamp, downHWrecamp];
%         HWpeakLat = [upHWpeakLat, downHWpeakLat]; % unused; HWLatency created and used from same variables later
        HWupslope = [upHWupslope, downHWupslope];
        HWdownslope = [upHWdownslope, downHWdownslope];
        
        %% mark events on the EEG
        disp('Marking events on the EEG dataset...')
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
                        eventsUp(nevt).duration = upHWperiod(n)*EEG.srate;
                    else
                        eventsUp(nevt).latency = latency + upHWpeakLat(n);
                        eventsUp(nevt).duration = 1; % one sample peak marker
                    end
                    eventsUp(nevt).channel = find(strcmp(ChName, ChOI{nch}));
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
                        eventsDown(nevt).duration = downHWperiod(n)*EEG.srate;
                    else
                        eventsDown(nevt).latency = latency + downHWpeakLat(n);
                        eventsDown(nevt).duration = 1; % one sample peak marker
                    end
                    eventsDown(nevt).channel = find(strcmp(ChName, ChOI{nch}));
                    eventsDown(nevt).urevent = [];
                    nevt = nevt + 1; % valid event counter
                end
                if n == length(segUp) % catch exceed end of array
                    break
                else
                    latency = latency + length(segUp{n+1}) + length(segDown{n});
                end
            end
        else
            latency = startDelay + length(segDown{1});
            nevt = 1; % valid event counter
            for n = 1:length(segUpDetect)
                if segUpDetect(n) == 1
                    eventsUp(nevt).type = eventName{1};
                    if isempty(peaks)
                        eventsUp(nevt).latency = latency;
                        eventsUp(nevt).duration = upHWperiod(n)*EEG.srate;
                    else
                        eventsUp(nevt).latency = latency + upHWpeakLat(n);
                        eventsUp(nevt).duration = 1; % one sample peak marker
                    end
                    eventsUp(nevt).channel = find(strcmp(ChName, ChOI{nch}));
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
                        eventsDown(nevt).duration = downHWperiod(n)*EEG.srate;
                    else
                        eventsDown(nevt).latency = latency + downHWpeakLat(n);
                        eventsDown(nevt).duration = 1; % one sample peak marker
                    end
                    eventsDown(nevt).channel = find(strcmp(ChName, ChOI{nch}));
                    eventsDown(nevt).urevent = [];
                    nevt = nevt + 1; % valid event counter
                end
                latency = latency + length(segUp{n}) + length(segDown{n});
            end
        end
        clear latency startDelay
        
        % concatenate event structures and sort
        EEG.event = [EEG.event eventsUp eventsDown];
        EEG = eeg_checkset(EEG, 'eventconsistency');
        EEG = eeg_checkset(EEG);
        
        %% label SW events with sleep stage
        disp('Finding sleep stage label for each SW event...')
        Event = EEG.event;
        
        if ~isfield(Event,'SleepStage')
            for i = 1:length(Event)
                Event(i).SleepStage = ''; % initialization
            end
        end
        
        evtIdx = find(ismember({Event.type},eventName));
        
        for iEvt = evtIdx % loop on event
            lastScoring = find(ismember({Event(1:iEvt).type},allSleepStages),1,'last');
            if ~isempty(lastScoring)
                Event(iEvt).SleepStage = Event(lastScoring).type;
            else
                Event(iEvt).SleepStage = '';
            end
        end
        EEG.event = Event;
        clear evtIdx iEvt
        
        %% populate table with results
        disp('Creating results table...')
        upN = [];
        upID = [];
        upCH = [];
        downN = [];
        downID = [];
        downCH = [];
        upSleepStage = [];
        downSleepStage =[];
        upEvt = 1; % valid event counter
        downEvt = 1; % valid event counter
        for n=1:length(segUpDetect)
            upN{n} = nfile;
            upID{n} = EEG.setname;
            upCH{n} = ChOI{nch};
            upSleepStage{n} = NaN;
            if segUpDetect(n) == 1
                upHWlatency{n} = eventsUp(upEvt).latency;
                upEvt = upEvt + 1; % valid event counter
            else
                upHWlatency{n} = NaN;
            end
        end
        for n=1:length(segDownDetect)
            downN{n} = nfile;
            downID{n} = EEG.setname;
            downCH{n} = ChOI{nch};
            downSleepStage{n} = NaN;
            if segDownDetect(n) == 1
                downHWlatency{n} = eventsDown(downEvt).latency;
                downEvt = downEvt + 1; % valid event counter
            else
                downHWlatency{n} = NaN;
            end
        end
        
        % concatenate
        HWlatency = [upHWlatency, downHWlatency];
        N = [upN, downN];
        ID = [upID, downID];
        CH = [upCH, downCH];
        sleepStage = [upSleepStage, downSleepStage];
        SWnew = table(N',ID',CH',sleepStage',HWlatency',HWperiod',HWfreq',HWpeak',HWintamp',HWrecamp',HWupslope',HWdownslope',...
            'VariableNames',{'N','ID','Channel','sleepStage','Latency','Duration','Frequency','PeakAmp','Area','AvgAmp','UpSlope','DownSlope'});
        SW = [SW ; SWnew]; % concatinate data for "all subjects" table
        
        % sort table by latency
        SW = sortrows(SW,[1 5]);
        
        clear upHWlatency downHWlatency HWlatency
        clear sleepStage segUpDetect segDownDetect segUp segDown idx_up idx_down
        clear upN downN upID downID upCH downCH upSleepStage downSleepStage N ID CH sleepStage HWlatency HWperiod HWfreq HWpeak HWintamp HWrecamp HWupslope HWdownslope HWpeakLat
        clear eventsUp eventsDown i n upEvt downEvt SWnew
        
    end
    
    %% mark for delete SW during movement artifact and outside NREM and outside Lights ON / Lights OFF tags
    disp('Deleting false positives...')
    Event = EEG.event;
    evtIdx = find(ismember({Event.type},eventName));
    badEvt = [];
    ToRmvArt = [];
    ToRmvSS = [];
    ToRmvLightsOff = [];
    ToRmvLightsOn = [];
    tagFoundFlag = 0;
    
    for iEvt = evtIdx % loop on event
        bad = find(ismember({Event(1:iEvt).type},badData),1,'last');
        if ~isempty(bad)
            if Event(bad).latency + Event(bad).duration > Event(iEvt).latency
                ToRmvArt(end+1) = iEvt;
            end
        end
    end
    for iEvt = evtIdx % loop on event
        bad = find(ismember({Event(1:iEvt).type},badSleepstages),1,'last');
        if ~isempty(bad)
            if Event(bad).latency + Event(bad).duration > Event(iEvt).latency
                ToRmvSS(end+1) = iEvt;
            end
        end
    end
    for iEvt = evtIdx % loop on event
        bad = find(ismember({Event(1:iEvt).type},lightsTags),1,'last');
        badEvt(end+1) = iEvt;
        if ~isempty(bad)
            if strcmp(Event(bad).type, lightsTags{1}) && tagFoundFlag == 0 % lights off
                    ToRmvLightsOff = badEvt; % to remove from first SW event up until lights off tag
                    tagFoundFlag = 1;
            elseif strcmp(Event(bad).type, lightsTags{2}) % lights on
                if Event(iEvt).latency > Event(bad).latency
                    ToRmvLightsOn(end+1) = iEvt; % to remove SW events after lights on tag
                end
            end
        end
    end
    clear tagFoundFlag bad badEvt Event evtIdx iEvt
    ToRmv = [ToRmvArt ToRmvSS ToRmvLightsOff ToRmvLightsOn];
    ToRmv = sort(ToRmv);
    
    % remove NaN values from the results table
    SW=SW(~any(ismissing(SW,NaN),2),:);
    
    % find the SW marked as bad in the table by their latency
    iRow=[];
    for nEvt=ToRmv
        row = find([SW.Latency{:}]' == EEG.event(nEvt).latency);
        if ~isempty(row)
            for i=1:length(row) % in case there are events with identical latencies
                iRow(end+1) = row(i);
            end
        end
    end
    iRow = unique(iRow); % to remove the repeat offenders
    
    % remove the bad events from the results table
    SW(iRow,:) = [];
    
    % now, lets remove the bad events from EEG.event
    EEG.event(ToRmv) = []; % delete
    
    % check that everything is ok
    EEG = eeg_checkset(EEG, 'eventconsistency');
    eeg_checkset(EEG);
    
    clear ToRmv ToRmvArt ToRmvSS ToRmvLightsOff ToRmvLightsOn
    
    %% find sleep stages for SW events and add to results table
    disp('Adding sleep stages to results table...')
    Event = EEG.event;
    evtIdx = find(ismember({Event.type},eventName));
    evtLatency = {Event(evtIdx).latency}';
    evtStage = {Event(evtIdx).SleepStage}';
    nevt = 1;
    for nRow = height(SW)-length(evtIdx)+1:height(SW)
        if contains(EEG.setname,char(table2cell(SW(nRow,2)))) && SW.Latency{nRow} == evtLatency{nevt}
            stage(nevt) = evtStage(nevt);
            nevt = nevt + 1;
        else
            error('Mismtach between table and SW events in EEG.event. Cannot assign sleep stage to SW.')
        end
    end
    SW(height(SW)-length(stage)+1:end,4) = stage'; % this is a little tricky, but it should do the job
    clear Event evtIdx evtLatency evtStage stage
    
    %% save out "single subject" table
    subidx = cell2mat(cellfun(@(n) strcmpi(n, EEG.setname), SW.ID, 'UniformOutput', 0)); % logical index of current subj's data
    SWsub = SW(subidx,:); 
    writetable(SWsub, [resultDir filesep EEG.setname '_SWevents.csv'], 'Delimiter', ',');
    
    %% save out the EEG dataset
    EEG.setname = [EEG.setname '_PAA'];
    disp(['Saving file ' EEG.setname '.set...'])
    pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', [resultDir filesep], 'savemode', 'onefile'); % for set files
    clear SWsub subidx SWnew startDelay nch datafilt ToRmv Event
    
end

%% save out the results
writetable(SW,[resultDir filesep 'AllSubjects_SWevents.csv'],'Delimiter',',');
save([resultDir filesep 'SWevents.mat'], 'SW');
clear nfile SW

disp('ALL DONE!!!')

toc

end