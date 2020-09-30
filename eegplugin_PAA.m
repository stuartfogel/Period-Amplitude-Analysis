function vers = eegplugin_PAA(fig, trystrs, catchstrs)

% eegplugin_detect_spindles() - EEGLAB plugin for period amplitude analysis
%
% Usage:
%   >> eegplugin_PAA(fig, trystrs, catchstrs)
%
% Inputs:
%   fig        - [integer] EEGLAB figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks. 
%
% See also:
%   PAA, pop_PAA, eeglab
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
%
% Copyright, Sleep Well. https://www.sleepwellpsg.com
%

vers = '1.6';
if nargin < 3
    error('eegplugin_PAA requires 3 arguments');
end

% add plugin folder to path
% -----------------------
if exist('pop_PAA.m','file')
    p = which('eegplugin_PAA');
    p = p(1:findstr(p,'eegplugin_PAA.m')-1);
    addpath(p);
end

% find tools menu
% ---------------------
menu = findobj(fig, 'tag', 'tools');

% menu callbacks
% --------------
PAA_cback = [ trystrs.no_check '[EEG,LASTCOM] = pop_PAA(EEG);' catchstrs.add_to_hist trystrs.no_check '[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);eeglab redraw;' catchstrs.add_to_hist ];

% create menus if necessary
% -------------------------
submenu = uimenu( menu, 'Label', 'Period Amplitude Analysis');
uimenu( submenu, 'Label', 'Period Amplitude Analysis', 'CallBack', PAA_cback);