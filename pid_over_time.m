%% PID_over_time_analysis.m
% Example script to compute Partial Information Decomposition (PID) over time
% using the MINT toolbox, ensuring the last dimension is "trials."

clearvars -except data_master
close all

% -------------------------------------------------------------------------
% Load data_master if it doesn't exist in workspace
% -------------------------------------------------------------------------
if ~exist('data_master', 'var')
    load('RCT_master_dataset_both_monkeys.mat'); % Adjust filename if needed
end
fprintf('Starting PID over time analysis...\n');

% -------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------
alignment_event = 'Align_to_cat_stim_on';
sample_rate     = 40000;     % 40 kHz
start_time      = 0.0;       % seconds
end_time        = 0.3;       % seconds
window_size     = 0.050;     % 50 ms window
window_step     = 0.010;     % 10 ms step
time_centers    = start_time : window_step : end_time;
n_windows       = length(time_centers);

% Subdivide each 50-ms window into smaller bins for MINT. E.g. 5 ms each:
sub_bin_size    = 0.005;     % 5 ms

% -------------------------------------------------------------------------
% Preprocess trial info and neurons
% -------------------------------------------------------------------------
session_ids     = [data_master.Bhv.session_id];
area_LIP        = 'MLIP';
area_FEF        = 'MFEF';
selected_date   = 20201211;  % Example session date
trial_info      = data_master.Bhv(selected_date == session_ids).Trial_info;

% Filter trial info (only correct, etc.)
params             = struct();
params.alignment   = alignment_event;
params.correct_only= 1;
[vert_trials, horizon_trials] = preprocess_trial_info(trial_info, params);

% Extract neurons for LIP, FEF
neur_LIP = data_master.Neuro.(area_LIP)( ...
             contains({data_master.Neuro.(area_LIP).NeuronID}, num2str(selected_date)));
neur_FEF = data_master.Neuro.(area_FEF)( ...
             contains({data_master.Neuro.(area_FEF).NeuronID}, num2str(selected_date)));

% -------------------------------------------------------------------------
% Initialize arrays to store PID results
% We'll compute one value per window (by averaging across sub-bins).
% For 2 sources (LIP, FEF) + 1 target (S), the 'PID_atoms' are:
% [Syn, Red, Unq1, Unq2]. We'll also get 'Joint' info as well.
% 
PID_Joint = nan(n_windows, 1);
PID_Syn   = nan(n_windows, 1);
PID_Red   = nan(n_windows, 1);
PID_Unq1  = nan(n_windows, 1);
PID_Unq2  = nan(n_windows, 1);

% -------------------------------------------------------------------------
% Loop over time windows
% -------------------------------------------------------------------------
for idx = 1:n_windows

    t_center = time_centers(idx);
    t_start  = max(t_center - window_size/2, 0);
    t_end    = t_center + window_size/2;
    
    fprintf('Processing window %d/%d: center=%.1f ms, [%.1f, %.1f] ms\n', ...
        idx, n_windows, t_center*1000, t_start*1000, t_end*1000);

    % 1) Extract spiking data in [t_start, t_end], subdivided by sub_bin_size:
    [S, X_LIP, Y_FEF] = extract_spike_matrices( ...
        vert_trials, neur_LIP, neur_FEF, sample_rate, ...
        [t_start, t_end], alignment_event, sub_bin_size);

    % 2) Force both arrays to have the same # of sub-bins (if needed)
    min_bins = min(size(X_LIP, 2), size(Y_FEF, 2));
    X_LIP = X_LIP(:, 1:min_bins, :);
    Y_FEF = Y_FEF(:, 1:min_bins, :);

    % 3) Permute to MINT's format *with trials in the LAST dimension*:
    %    The original shape: [nTrials, nTimeBins, nNeurons]
    %    MINT wants:         [nNeurons, nTimeBins, nTrials]
    X1 = permute(X_LIP, [3, 2, 1]);
    X2 = permute(Y_FEF, [3, 2, 1]);

    % 4) Reshape S so the last dimension is #trials (e.g. [1 x 1 x nTrials])
    nTrials = length(S);
    S = reshape(S, [1, 1, nTrials]);

    fprintf('  Size X1: %s | Size X2: %s | #Trials: %d\n', ...
        mat2str(size(X1)), mat2str(size(X2)), nTrials);

    % Basic checks
    if size(X1,2) ~= size(X2,2)
        fprintf('  Skipping: Mismatch in # time bins!\n');
        continue
    end
    if size(X1,3) ~= size(X2,3) || size(X1,3) ~= size(S,3)
        fprintf('  Skipping: Mismatch in # trials!\n');
        continue
    end
    if size(X1,2) < 1
        fprintf('  Skipping: Not enough sub-bins!\n');
        continue
    end
    if size(X1,3) < 30
        fprintf('  Skipping: Not enough trials!\n');
        continue
    end

    % 5) PID options
    nTimeBins = size(X1,2);
    PID_opts = struct();
    PID_opts.tpres              = {nTimeBins};  % 'present' time for each source
    PID_opts.tau                = {min(2, nTimeBins-1)}; % small tau
    PID_opts.redundancy_measure = 'I_min';      % e.g. 'I_BROJA', 'I_min', 'I_MMI'
    PID_opts.bin_method         = {'eqpop','eqpop','none'};
    PID_opts.n_bins             = {3, 3};
    PID_opts.bias               = 'shuffSub';   % or 'plugin', 'infoCorr', etc.
    PID_opts.pid_constrained    = true;
    PID_opts.supressWarnings    = true;
    PID_opts.computeNulldist    = false;
    PID_opts.parallel_sampling  = true;

    % If tpres <= tau or # timebins too small, skip
    if PID_opts.tpres{1} <= PID_opts.tau{1}
        fprintf('  Skipping: tpres <= tau.\n');
        continue
    end

    % 6) Call PID
    try
        inputs     = {X1, X2, S};  
        % We request 'Joint' (I(X1,X2;S)) and 'PID_atoms' ([Syn, Red, Unq1, Unq2]):
        outputList = {'Joint','PID_atoms'};
        [PID_corrected, ~, ~] = PID(inputs, outputList, PID_opts);

        % 7) MINT returns a value for each sub-bin => e.g. 1 x nTimeBins
        % For 'PID_atoms', we get a 4 x nTimeBins array (rows: Syn,Red,Unq1,Unq2).
        % Typically we average across sub-bins to get a single measure for that window:
        joint_vals = PID_corrected{1};   % shape [1 x nTimeBins]
        atom_vals  = PID_corrected{2};   % shape [4 x nTimeBins]

        % Average across sub-bins:
        Joint      = mean(joint_vals, 2);     % => scalar
        Syn        = mean(atom_vals(1,:), 2);
        Red        = mean(atom_vals(2,:), 2);
        Unq1       = mean(atom_vals(3,:), 2);
        Unq2       = mean(atom_vals(4,:), 2);

        % Store
        PID_Joint(idx) = Joint;
        PID_Syn(idx)   = Syn;
        PID_Red(idx)   = Red;
        PID_Unq1(idx)  = Unq1;
        PID_Unq2(idx)  = Unq2;

        fprintf('  PID Joint=%.4f, Syn=%.4f, Red=%.4f, Unq1=%.4f, Unq2=%.4f\n', ...
            Joint, Syn, Red, Unq1, Unq2);

    catch ME
        fprintf('  Error during PID computation: %s\n', ME.message);
        continue
    end
end

fprintf('Finished PID over time analysis.\n');

% -------------------------------------------------------------------------
% Plot the results
% -------------------------------------------------------------------------
% We'll plot Joint, Syn, Red, Unq1, Unq2 over the n_windows
figure('Name','PID Over Time','Color','w');
hold on;
plot(time_centers*1000, PID_Joint, 'k-', 'LineWidth', 2, 'DisplayName','Joint');
plot(time_centers*1000, PID_Syn,   'r-', 'LineWidth', 2, 'DisplayName','Syn');
plot(time_centers*1000, PID_Red,   'b-', 'LineWidth', 2, 'DisplayName','Red');
plot(time_centers*1000, PID_Unq1,  'g-', 'LineWidth', 2, 'DisplayName','Unq1');
plot(time_centers*1000, PID_Unq2,  'm-', 'LineWidth', 2, 'DisplayName','Unq2');
xlabel('Time from Stimulus Onset (ms)');
ylabel('Partial Information (bits)');
title('PID between LIP and FEF Over Time');
legend('Location','best');
grid on; box off;
ylim([-0.005, 0.01]); % adjust as needed