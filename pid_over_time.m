%% Load data
clearvars -except data_master
close all

if ~exist('data_master', 'var')
    load('RCT_master_dataset_both_monkeys.mat');
end

fprintf('Starting FIT over time analysis...\n');

%% Settings
alignment_event = 'Align_to_cat_stim_on';
sample_rate = 40000; % 40kHz
start_time = 0.0;    % seconds
end_time = 0.3;      % seconds
window_size = 0.050; % 50 ms window
window_step = 0.010; % 10 ms step
time_centers = start_time:window_step:end_time;
n_windows = length(time_centers);

% Initialize FIT results
FIT_AtoB = nan(n_windows, 1);
FIT_BtoA = nan(n_windows, 1);

%% Preprocess trial info and neurons
session_ids = [data_master.Bhv.session_id];
area_LIP = 'MLIP'; 
area_FEF = 'MFEF'; 
selected_date = 20201211; % Example session

trial_info = data_master.Bhv(selected_date == session_ids).Trial_info;

% Preprocess trial info
params = struct();
params.alignment = alignment_event;
params.correct_only = 1;

[vert_trials, horizon_trials] = preprocess_trial_info(trial_info, params);

% Extract neurons
neur_LIP = data_master.Neuro.(area_LIP)(contains({data_master.Neuro.(area_LIP).NeuronID}, num2str(selected_date)));
neur_FEF = data_master.Neuro.(area_FEF)(contains({data_master.Neuro.(area_FEF).NeuronID}, num2str(selected_date)));

%% Loop over time windows
for idx = 1:n_windows

    t_center = time_centers(idx);
    t_start = max(t_center - window_size/2, 0);
    t_end = t_center + window_size/2;
    
    fprintf('Processing window %d/%d: center=%.1f ms, [%.1f, %.1f] ms\n', ...
        idx, n_windows, t_center*1000, t_start*1000, t_end*1000);

    % Extract spike counts
    [S, X_LIP, Y_FEF] = extract_spike_matrices(vert_trials, neur_LIP, neur_FEF, sample_rate, [t_start, t_end], alignment_event);

    % 🛠 NEW: Force-align number of time bins
    min_bins = min(size(X_LIP,2), size(Y_FEF,2));
    X_LIP = X_LIP(:,1:min_bins,:);
    Y_FEF = Y_FEF(:,1:min_bins,:);

    % Reshape to MINT format
    X1 = permute(X_LIP, [3 2 1]); % [trials, time, neurons]
    X2 = permute(Y_FEF, [3 2 1]);
    S = S(:)'; % Labels into row vector

    % Debug prints
    fprintf('  Size X1: %s | Size X2: %s | #Trials: %d\n', mat2str(size(X1)), mat2str(size(X2)), length(S));

    % Check if data is valid
    if size(X1,2) <= 2
        fprintf('  Skipping: Not enough time bins (only %d bins).\n', size(X1,2));
        continue
    end
    if size(X1,3) <= 30
        fprintf('  Skipping: Not enough trials (only %d trials).\n', size(X1,3));
        continue
    end
    if size(X2,2) ~= size(X1,2) || size(X2,3) ~= size(X1,3)
        fprintf('  Skipping: Mismatch in X1 and X2 shapes!\n');
        continue
    end

    % Setup FIT options
    window_bins = size(X1,2);
    FIT_opts = struct();
    FIT_opts.tpres = {window_bins};
    FIT_opts.tau = {min(5, window_bins-2)};
    if FIT_opts.tpres{1} <= FIT_opts.tau{1}
        fprintf('  Skipping: tpres <= tau.\n');
        continue
    end
    FIT_opts.redundancy_measure = 'I_min';
    FIT_opts.bin_method = {'eqpop', 'eqpop', 'none'};
    FIT_opts.n_bins = {3, 3};
    FIT_opts.bias = 'shuffSub';
    FIT_opts.pid_constrained = true;
    FIT_opts.supressWarnings = true;
    FIT_opts.computeNulldist = false;
    FIT_opts.parallel_sampling = true;

    % Try to compute FIT
    try
        inputs = {X1, X2, S};
        outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'};
        [FIT_corrected, ~, ~] = FIT(inputs, outputList, FIT_opts);

        % Store results
        FIT_AtoB(idx) = FIT_corrected{1};
        FIT_BtoA(idx) = FIT_corrected{2};
        fprintf('  FIT A->B: %.5f | FIT B->A: %.5f\n', FIT_AtoB(idx), FIT_BtoA(idx));

    catch ME
        fprintf('  Error during FIT computation: %s\n', ME.message);
        continue
    end
end

fprintf('Finished FIT over time analysis.\n');

%% Plot the results
figure;
hold on;
plot(time_centers*1000, FIT_AtoB, 'r-', 'LineWidth', 2); % LIP -> FEF
plot(time_centers*1000, FIT_BtoA, 'b-', 'LineWidth', 2); % FEF -> LIP
xlabel('Time from Stimulus Onset (ms)');
ylabel('Feature Information Transfer (bits)');
title('FIT between LIP and FEF Over Time');
legend('LIP \rightarrow FEF', 'FEF \rightarrow LIP');
ylim([-0.005, 0.005]); % Adjust based on real results
grid on;
box off;