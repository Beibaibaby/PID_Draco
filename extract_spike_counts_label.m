function [S, X, Y] = extract_spike_counts_label( ...
         trials, neur_LIP, neur_FEF, Fs, win, align_event, label_field)
% Count spikes in the window `win` (s) relative to `align_event`.
%
% RETURNS
%   S   – trial label vector (category / direction / etc.)
%   X   – [trials × nLIP]   spike counts
%   Y   – [trials × nFEF]   spike counts
%
% EXAMPLE
%   [S,X,Y] = extract_spike_counts_label(trials,neur_LIP,neur_FEF,40e3,...
%                                       [0 0.05],'Align_to_cat_stim_on',...
%                                       'category');

if nargin < 7, label_field = 'direction'; end

nT    = numel(trials);
nLIP  = numel(neur_LIP);
nFEF  = numel(neur_FEF);

X = zeros(nT,nLIP,'single');
Y = zeros(nT,nFEF,'single');
S = zeros(nT,1,'int8');

for k = 1:nT
    tr  = trials(k);
    t0  = tr.(align_event);              % alignment time in samples
    t0s = t0 / Fs;                       % convert to seconds

    S(k) = tr.(label_field);             % category label (-1 / 1)

    % LIP units
    for u = 1:nLIP
        spk = neur_LIP(u).NeuronSpkT / Fs;
        X(k,u) = sum( spk >= t0s + win(1) & spk < t0s + win(2) );
    end

    % FEF units
    for u = 1:nFEF
        spk = neur_FEF(u).NeuronSpkT / Fs;
        Y(k,u) = sum( spk >= t0s + win(1) & spk < t0s + win(2) );
    end
end
end