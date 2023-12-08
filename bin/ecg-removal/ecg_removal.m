function [cleaned_ats] = ecg_removal(signal)
addpath('filters'); addpath('ecg_utils'); addpath('template_subtraction');
use_filtfilt = true;

% signal should be in one row

fs = 1000; fpl = 50;  % fpl is powerline frequency (typically 50Hz or 60Hz)

% power line interference removal
signal = butter_filt_stabilized(signal, [fpl-1 fpl+1], fs, 'stop', use_filtfilt, 2);

% mild high-pass filtering (two versions, one for R peak detection and one for cardiac artifact removal) 
% to remove baseline wander and some ECG interference (especially in the 20Hz version)
signalhp05 = butter_filt_stabilized(signal, 5, fs, 'high', use_filtfilt, 6);
signalhp20 = butter_filt_stabilized(signal, 20, fs, 'high', use_filtfilt, 6);

% R peak detection, slightly modified version of Pan-Tompkins
rpeaks = peak_detection(signalhp05, fs);

% This is the actual cardiac artifact removal step
cleaned_ats = adaptive_template_subtraction(signalhp20, rpeaks, fs);
% Wavelet denoising is another very robust alternative
cleaned_swt = swtden(signal, rpeaks, fs, 'h', 3, 'db2', 4.5);
% Depending on the use case, even a simple HP100 might do
cleaned_hp100 = butter_filt_stabilized(signal, 100, fs, 'high', use_filtfilt, 6);
