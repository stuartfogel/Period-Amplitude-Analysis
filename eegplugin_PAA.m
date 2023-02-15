function vers = eegplugin_PAA(fig, trystrs, catchstrs)

% eegplugin_PAA() - EEGLAB plugin for period amplitude analysis
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
% May 26, 2022 Revised 1.9: corrected bug in calculation on integrated
%   amplitude - SF
% Aug 30, 2022 Revised 2.0: new SW parameters added to event structure. 
% Dec 18, 2022 Revised 2.01: better handling of GUI input with multiple channels. 
% Dec 28, 2022 Revised 2.02: modified menu layout
% Jan 21, 2023 Revised 2.1: batch uses GUI
% Feb 12, 2023 Revised 2.1.1: added progress bars
% Feb 15, 2023 Revised 2.1.2: handle multiple bad data markers
%
% Copyright, Sleep Well. https://www.sleepwellpsg.com
%

vers = '2.1.2';
if nargin < 3
    error('eegplugin_PAA requires 3 arguments');
end

% add plugin folder to path
if exist('pop_PAA.m','file')
    p = which('eegplugin_PAA');
    p = p(1:findstr(p,'eegplugin_PAA.m')-1);
    addpath(p);
end

% find tools menu
menu = findobj(fig, 'tag', 'tools');

% menu callbacks
PAA_cback = [ trystrs.no_check '[EEG,LASTCOM] = pop_PAA(EEG);' catchstrs.add_to_hist trystrs.no_check '[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);eeglab redraw;' catchstrs.add_to_hist ];

% create menu
uimenu(menu, 'Label', 'Period Amplitude Analysis', 'CallBack', PAA_cback);
