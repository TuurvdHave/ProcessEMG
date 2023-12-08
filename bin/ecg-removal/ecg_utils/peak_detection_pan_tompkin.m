function [peaks, polarity] = peak_detection_pan_tompkin(data, fs, varargin)
%PEAK_DETECTION_PAN_TOMPKIN R-peak detector based on Pan-Tompkins method.
%
% inputs:
% x: vector of input data
% fs: sampling rate
% wlen: moving average window length (default = 150ms)
% fp1: lower cut-off frequency (default = 10Hz)
% fp2: upper cut-off frequency (default = 33.3Hz)
% th: detection threshold (default = 0.2)
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% peaks: vector of R-peak impulse train
%
% Original version of this file:
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Original work Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% Modified Eike Petersen, Institute for Electrical Engineering in Medicine,
% University of Lübeck.
%
% Modifications:
% * Use a moving / adaptive threshold instead of a fixed threshold
% * More robust signal polarity detection
% * Use Matlab movmedian implementation instead of custom implementation
% * Used more robust peak detection in the last step

if(nargin > 2  && ~isempty(varargin{1})),
    winlen = varargin{1};
else
    winlen = .150; % 150ms
end

if(nargin > 3 && ~isempty(varargin{2})),
    fp1 = varargin{2};
else
    fp1 = 10;
end

if(nargin > 4  && ~isempty(varargin{3})),
    fp2 = varargin{3};
else
    fp2 = 33.3;
end

if(nargin > 5  && ~isempty(varargin{4})),
    thr = varargin{4};
else
    thr = 0.2;
end

if(nargin > 6  && ~isempty(varargin{5})),
    flag = varargin{5};
else
    % MODIFICATION BY IME
    % Determine polarity of signal: do we expect R peaks to be maxima or
    % minima?
    maxima = movmax(data, 2*fs);
    minima = movmin(data, 2*fs);
    flag = sum(maxima > abs(minima)) > length(data)/2;
end
polarity = flag;

if polarity
    peak_sign = 1;
else
    peak_sign = -1;
end

N = length(data);
data = data(:);

L1 = round(fs/fp2);    % first zero of the LP filter is placed at f = 33.3Hz;
L2 = round(fs/fp1);    % first zero of the HP filter is placed at f = 3Hz;

% MODIFICATION BY IME: use Matlab implementation of moving median
x0 = data - movmedian(data, round(fs*winlen/3));

% LP filter
x = filter([1 zeros(1,L1-1) -1],[L1 -L1],x0);
x = filter([1 zeros(1,L1-1) -1],[L1 -L1],x);
x = [x(L1:end);zeros(L1-1,1) + x(end)]; % lag compensation

% HP filter
y = filter([L2-1 -L2 zeros(1,L2-2) 1],[L2 -L2],x);

% differentiation
z = diff([y(1) ; y]);

% squaring
w = z.^2;

% moving average
L3 = round(fs*winlen);
v = filter([1 zeros(1,L3-1) -1],[L3 -L3],w);
v = [v(round(L3/2):end);zeros(round(L3/2)-1,1) + v(end)]; % lag compensation

% MODIFICATION BY IME
% New, moving threshold calculation
% Calculate max over moving ~1 beat window
vmovmax = movmax(v, round(fs*1.5), 'omitnan');
% Calculate mean over moving 3s window, corresponding to the average maxima
% of maybe the neighboring 5 beats
vmovmaxmean = movmean(vmovmax, round(fs*3));
p = v > (thr*vmovmaxmean);

% find peaks in v
[~, vpeakloc] = findpeaks(v, 'MinPeakDistance', round(winlen*fs*2));
% Reject peaks below threshold
vpeakloc = vpeakloc(p(vpeakloc));

% Now find corresponding peaks in the raw signal
% edge detection
rising  = find(diff([0 ; p])==1);      % rising edges
falling = find(diff([p ; 0])==-1);     % falling edges

if( length(rising) == length(falling)-1 )
    rising = [1 ; rising];
elseif( length(rising) == length(falling)+1 )
    falling = [falling ; N];
end

peakloc = zeros(length(vpeakloc), 1);

for i=1:length(vpeakloc)
	rising_edge = max(rising(rising < vpeakloc(i)));
	falling_edge = min(falling(falling > vpeakloc(i)));
	[~, mx] = max(peak_sign * data(rising_edge:falling_edge) );
	peakloc(i) = mx - 1 + rising_edge;
end

if peakloc(1) == 1
    peakloc = peakloc(2:end);
end

if peakloc(end) == length(data)
    peakloc = peakloc(1:end-1);
end

peaks = zeros(1,N);
peaks(peakloc) = 1;