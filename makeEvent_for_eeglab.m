function event = makeEvent_for_eeglab(sessionMetaBehav)
% This code convert sessionMetaBehav table to triggers, each trial will
% create Init and Stim triggers, with trial type labeled eg. 'StimCorD4'

%%
addpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\toolboxes\eeglab2019_0');

% type latency urevent
colDD = 6;
colHit = 9;
colInit = 16;
colStim = 17;
colTouch = 18;

nTrial = size(sessionMetaBehav,1);
correctTypes = {'Inc', 'Cor', 'Pre', 'Omi'};
alignTypes   = {'Init', 'Stim'};
delayDurations = [4,5,6];

for iTrial = 1 : nTrial
    delayType = num2str(sessionMetaBehav(iTrial,colDD).DelayDuration);
    correctType = correctTypes{sessionMetaBehav(iTrial,colHit).HitMiss+1};
    event(iTrial).type = ['Init' correctType 'D' delayType];
    event(iTrial+nTrial).type = ['Stim' correctType 'D' delayType];
    
    event(iTrial).latency=sessionMetaBehav(iTrial,colInit).Init*1000; % to msec
    event(iTrial+nTrial).latency=sessionMetaBehav(iTrial,colStim).StimOnset*1000; % to msec
    
    event(iTrial).urevent = iTrial; 
    event(iTrial+nTrial).urevent = iTrial;    
end
