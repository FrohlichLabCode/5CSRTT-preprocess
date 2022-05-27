%
%
% LFP data pre-processing
%
% Sangtae Ahn (sangtae_ahn@med.unc.edu)
% Frohlich Lab.
%
% first written by Sangtae Ahn on July 18, 2019
%
%

close all;
clear;
clc;



%% 
addpath('eeglab2019_0');
eeglab;
pop_editoptions( 'option_savetwofiles', 1,'option_single', 0);

load('lfpMat.mat');

highCut = 300;

%%  
EEG.data = lfpMat(:,1:500000);
EEG.srate = lfpFs;
EEG = eeg_checkset(EEG);
eeglab redraw
pop_saveset(EEG,'filename',['lfpMat.set']);

% %% spectra
% pop_eegplot(EEG,1,1,1); % plot raw traces
% 
% figure; 
% pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',[1 500],'electrodes','off');
% 
% EEG = pop_eegfiltnew(EEG, 'hicutoff',300);
% figure; 
% pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',[1 500],'electrodes','off');

%% removal Line noise
addpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeSangtae\LFP\removeLineNoise_SpectrumEstimation');

LFP_wave = removeLineNoise_SpectrumEstimation(EEG.data, EEG.srate, ...
    'LF = 31, NH = 10, M = 1024, HW = 3');
%       NH  = number of harmonics. (default: 1)
%       LF  = line frequency. (default: 50 Hz)
%       TOL = error tolerence. (default: 1% tol)
%       HW  = half-width of peak in samples. (default: 2)
%       M   = Size of window. (For fs = 1 kHz, m = 1024). If M is set, TOL
%                              has no effect.
EEG.data=LFP_wave;
% figure;
% pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',[1 500],'electrodes','off');

[spec1 freq]=pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',[1 500],'electrodes','off','plot','off');
% figure;plot(freq,spec1);

ROI = ['PFC'; 'LPl'; 'PPC'; 'VIS'];
figure;
for iROI = 1 : length(ROI)
    subplot(2,2,iROI);
    plot(freq,spec1(iROI*16-15:iROI*16,:)');
    title(ROI(iROI,:));
    xlabel('frequency (Hz)');
    ylabel('power (dB)');
end






