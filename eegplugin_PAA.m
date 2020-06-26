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
% Author: 
%   Stuart Fogel <sfogel@uottawa.ca>
%   School of Psychology
%   University of Ottawa
%
% Copyright (C) <2020>  Stuart Fogel, http://socialsciences.uottawa.ca/sleep-lab/
%

vers = '1.0';
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