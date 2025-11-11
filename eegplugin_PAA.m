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
%
% Copyright, Sleep Well. https://www.sleepwellpsg.com
%

vers = '2.1.6'; % handles NaN in data and other minor fixes.
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
