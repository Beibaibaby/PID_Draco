%% Load data
clearvars -except data_master
close all

if ~exist('data_master','var')
    load('RCT_master_dataset_both_monkeys.mat');
end

fprintf('Starting CATEGORY-specific FIT over time analysis...\n');

%% Settings
alignment_event = 'Align_to_cat_stim_on';
sample_rate = 40000; % 40kHz
start_time = 0.0;    % seconds
end_time = 0.3;      % seconds
window_size = 0.030; % 50 ms window
window_step = 0.005; % slide every 10 ms
time_centers = start_time:window_step:end_time;
n_windows = length(time_centers);

% Initialize storage
FIT_AtoB = nan(n_windows, 1);
FIT_BtoA = nan(n_windows, 1);

%% Preprocess trial info
session_ids = [data_master.Bhv.session_id];
brain_areas = fieldnames(data_master.Neuro);
area_LIP = 'MLIP'; 
area_FEF = 'MFEF'; 

selected_date = 20201211; % Example session
trial_info = data_master.Bhv(selected_date == session_ids).Trial_info;

params = struct();
params.alignment = alignment_event;
params.correct_only = 1;

[vert_trial_info, horizon_trial_info] = preprocess_trial_info(trial_info, params);

neur_info_LIP = data_master.Neuro.(area_LIP)(contains({data_master.Neuro.(area_LIP).NeuronID}, num2str(selected_date)));
neur_info_FEF = data_master.Neuro.(area_FEF)(contains({data_master.Neuro.(area_FEF).NeuronID}, num2str(selected_date)));

trial_info_use = vert_trial_info; % Use vertical trials for now

%% Loop over time windows
for idx = 1:n_windows
    t_center = time_centers(idx);
    t_start = max(t_center - window_size/2, 0);
    t_end = t_center + window_size/2;

    fprintf('Processing window %d/%d: center=%.1f ms, [%.1f, %.1f] ms\n', ...
        idx, n_windows, t_center*1000, t_start*1000, t_end*1000);

    % Extract spike counts, and extract S_temp (very important)
    [S_temp, X_LIP, Y_FEF] = extract_spike_matrices(trial_info_use, neur_info_LIP, neur_info_FEF, sample_rate, [t_start, t_end], alignment_event);

    % Correct: use S_temp for categories
    S_cat = S_temp; % Use the stimulus features returned
    S_cat(S_cat == -1) = 2; % Remap -1 -> 2
    S_cat = S_cat(:)'; % Ensure row vector

    % Reshape to MINT input format
    X1 = permute(X_LIP, [3 2 1]); % [1, time, trials]
    X2 = permute(Y_FEF, [3 2 1]);

    % Safety check
    if size(X1,2) <= 2 || size(X1,3) <= 30
        fprintf('  Skipping window: too small.\n');
        continue
    end

    if length(S_cat) ~= size(X1,3)
        fprintf('  Warning: Length of S_cat (%d) does not match number of trials (%d)! Skipping.\n', length(S_cat), size(X1,3));
        continue
    end

    % Setup FIT
    inputs = {X1, X2, S_cat};
    outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'};
    window_bins = size(X1,2);

    FIT_opts = struct();
    FIT_opts.tpres = {window_bins};
    FIT_opts.tau = {min(3, window_bins-2)};
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
        [FIT_corrected, ~, ~] = FIT(inputs, outputList, FIT_opts);
        FIT_AtoB(idx) = FIT_corrected{1}; % LIP->FEF
        FIT_BtoA(idx) = FIT_corrected{2}; % FEF->LIP
        fprintf('  FIT A->B: %.5f | B->A: %.5f\n', FIT_AtoB(idx), FIT_BtoA(idx));
    catch ME
        fprintf('  Error at idx %d: %s\n', idx, ME.message);
        continue
    end
end

fprintf('Finished CATEGORY-specific FIT analysis.\n');

%% Plot the results
figure;
hold on;
plot(time_centers*1000, FIT_AtoB, 'r-', 'LineWidth', 2);
plot(time_centers*1000, FIT_BtoA, 'b-', 'LineWidth', 2);
xlabel('Time from Stimulus Onset (ms)');
ylabel('Category Information Transfer (bits)');
title('Category-specific FIT between LIP and FEF Over Time');
legend('LIP \rightarrow FEF', 'FEF \rightarrow LIP');
ylim([-0.005, 0.005]);
grid on;
box off;