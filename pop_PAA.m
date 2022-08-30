function [EEG,com] = pop_PAA(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EEGlab GUI for Period Amplitude Analysis (PAA) on EEGlab data to
% detect and measure slow waves. 
%
% Copyright, Sleep Well. https://www.sleepwellpsg.com
%
% Method adpated from:
%
% Feinberg et al, 1978. https://doi.org/10.1016/0013-4694(78)90266-3
% Geering et al, 1993. https://doi.org/10.1111/j.1365-2869.1993.tb00074.x
% Bersagliere and Achermann, 2010. https://doi.org/10.1111/j.1365-2869.2009.00775.x
% ...also, see recent (2020's-ish) papers by Tononi's group
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