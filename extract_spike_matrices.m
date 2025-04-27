function [S, X, Y] = extract_spike_matrices(trials, neur_LIP, neur_FEF, sample_rate, window, alignment)

n_trials = length(trials);

% Preallocate
X = zeros(n_trials, length(neur_LIP));
Y = zeros(n_trials, length(neur_FEF));
S = zeros(n_trials, 1);

for i = 1:n_trials
    trial = trials(i);
    align_time = trial.(alignment);

    % Direction label
    S(i) = trial.direction; 

    % LIP
    for j = 1:length(neur_LIP)
        spk = neur_LIP(j).NeuronSpkT / sample_rate;
        X(i,j) = sum(spk >= (align_time + window(1)) & spk <= (align_time + window(2)));
    end

    % FEF
    for j = 1:length(neur_FEF)
        spk = neur_FEF(j).NeuronSpkT / sample_rate;
        Y(i,j) = sum(spk >= (align_time + window(1)) & spk <= (align_time + window(2)));
    end
end
end