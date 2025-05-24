function [S, X, Y] = extract_spike_matrices(trials, neur_LIP, neur_FEF, ...
    sample_rate, window, alignment, bin_size)
% extract_spike_matrices  Returns binned spike counts for LIP and FEF 
%   in the specified time window [window(1), window(2)], subdivided into 
%   smaller bins of 'bin_size' length.
%
%   X: [nTrials, nBins, nNeuronsLIP]
%   Y: [nTrials, nBins, nNeuronsFEF]
%   S: [nTrials, 1] (labels; e.g. direction)

if nargin < 7
    % fallback if the user doesn't provide a sub-bin size
    bin_size = 0.005; % 5 ms default
end

n_trials  = length(trials);
nNeurLIP  = length(neur_LIP);
nNeurFEF  = length(neur_FEF);
S         = zeros(n_trials, 1);

total_dur = window(2) - window(1);
n_bins    = floor(total_dur / bin_size);

% Pre-allocate X, Y with third dimension = # of neurons
X = zeros(n_trials, n_bins, nNeurLIP);
Y = zeros(n_trials, n_bins, nNeurFEF);

for i = 1:n_trials
    trial      = trials(i);
    align_time = trial.(alignment);

    % Direction or other label
    S(i) = trial.direction;

    % For each LIP neuron
    for j = 1:nNeurLIP
        spk_times = double(neur_LIP(j).NeuronSpkT) / sample_rate;

        % Loop over sub-bins
        for b = 1:n_bins
            t0 = align_time + window(1) + (b-1)*bin_size;
            t1 = t0 + bin_size;
            X(i,b,j) = sum(spk_times >= t0 & spk_times < t1);
        end
    end

    % For each FEF neuron
    for j = 1:nNeurFEF
        spk_times = double(neur_FEF(j).NeuronSpkT) / sample_rate;

        % Loop over sub-bins
        for b = 1:n_bins
            t0 = align_time + window(1) + (b-1)*bin_size;
            t1 = t0 + bin_size;
            Y(i,b,j) = sum(spk_times >= t0 & spk_times < t1);
        end
    end
end

end