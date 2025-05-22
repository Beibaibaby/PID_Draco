function fit_dir_cat_allSessions_a()
% FIT_DIR_CAT_ALLSESSIONS_CORRECTED
%
% This script loops over all sessions, computes FIT for direction and
% category labels within time windows for each session.
% **Modifications:**
%   - Saves results (.mat) and plots (.png) for EACH session individually
%     after processing is complete for that session.
%   - Checks if a session's .mat file already exists and skips computation,
%     loading the previous results instead (simple resume capability).
%   - Averages the results across all processed/loaded sessions at the end.
%   - Uses correct data shaping and MINT options.

%% 0) Setup Output Directories
output_mat_dir = 'session_results_mat';
output_fig_dir = 'session_results_figs';
if ~exist(output_mat_dir, 'dir')
   mkdir(output_mat_dir);
   fprintf('[INFO] Created directory for session .mat files: %s\n', output_mat_dir);
end
if ~exist(output_fig_dir, 'dir')
   mkdir(output_fig_dir);
   fprintf('[INFO] Created directory for session figures: %s\n', output_fig_dir);
end

%% 1) Ensure data_master is loaded
if ~exist('data_master','var')
    try
        load('RCT_master_dataset_both_monkeys.mat');
        fprintf('[INFO] data_master loaded from file.\n');
    catch ME
        error('Failed to load data_master: %s', ME.message);
    end
else
    fprintf('[INFO] data_master already in workspace.\n');
end

%% 2) Identify all sessions
if ~isfield(data_master, 'Bhv') || ~isfield(data_master.Bhv, 'session_id')
    error('data_master does not contain Bhv.session_id field.');
end
allSessionIDs = [data_master.Bhv.session_id];
uniqueSessions = unique(allSessionIDs);
nSessions = length(uniqueSessions);
fprintf('[INFO] Found %d unique sessions.\n', nSessions);

%% 3) Parse NeuronIDs -> session dates for LIP & FEF
area_LIP = 'MLIP';
area_FEF = 'MFEF';
if ~isfield(data_master, 'Neuro') || ~isfield(data_master.Neuro, area_LIP) || ~isfield(data_master.Neuro, area_FEF)
    error('data_master.Neuro missing LIP (%s) or FEF (%s) data.', area_LIP, area_FEF);
end
allNeurLIP = data_master.Neuro.(area_LIP);
allNeurFEF = data_master.Neuro.(area_FEF);
NeurIDs_LIP = {allNeurLIP.NeuronID};
NeurIDs_FEF = {allNeurFEF.NeuronID};
sess_dates_LIP = parseSessionDates(NeurIDs_LIP);
sess_dates_FEF = parseSessionDates(NeurIDs_FEF);

%% 4) MINT parameters and time window settings
alignment_event = 'Align_to_cat_stim_on';
sample_rate     = 40000;  % 40 kHz
window_size     = 0.050;  % 50 ms window (Adjusted back from user example)
window_step     = 0.010;  % 20 ms step (Adjusted back from user example)
start_time      = 0.0;    % seconds relative to alignment event
end_time        = 0.3;    % seconds relative to alignment event
time_centers    = start_time : window_step : end_time;
n_windows       = length(time_centers);
sub_bin_size    = 0.005;  % 2 ms sub-bins for MINT analysis
min_trials_thresh = 30;   % Minimum trials needed per session/window

% --- FIT Default Options ---
FIT_opts_base = struct();
FIT_opts_base.tauVal              = 2; % 2 bins * 5ms/bin = 10ms delay
FIT_opts_base.redundancy_measure = 'I_min';
FIT_opts_base.bin_method         = {'eqpop','eqpop','none'};
FIT_opts_base.n_bins             = {3,3};
FIT_opts_base.bias               = 'plugin';
FIT_opts_base.xtrp               = 10;
FIT_opts_base.shuff              = 20;
FIT_opts_base.pid_constrained    = true;
FIT_opts_base.supressWarnings    = true;
FIT_opts_base.computeNulldist    = false;
FIT_opts_base.n_samples          = 100;
FIT_opts_base.parallel_sampling  = false;
FIT_opts_base.shuffling          = {'A', 'B'};
% --------------------------

% We'll store results across sessions
all_dirAtoB = nan(n_windows, nSessions);
all_dirBtoA = nan(n_windows, nSessions);
all_catAtoB = nan(n_windows, nSessions);
all_catBtoA = nan(n_windows, nSessions);
sessions_processed_count = 0; % Counter for successfully processed/loaded sessions

%% 5) Loop over sessions
for sIdx = 1:nSessions
    thisDate = uniqueSessions(sIdx);
    fprintf('\n=== Session %d (%d/%d) ===\n', thisDate, sIdx, nSessions);

    % --- Check for existing results file ---
    session_mat_filename = fullfile(output_mat_dir, sprintf('session_%d_results.mat', thisDate));
    if exist(session_mat_filename, 'file')
        fprintf('  [INFO] Results file found: %s. Loading previous results.\n', session_mat_filename);
        try
            load(session_mat_filename, 'sessionResult');
            % Check if loaded data structure is compatible
            if isfield(sessionResult, 'dirAtoB') && isfield(sessionResult,'time_centers') && ...
               length(sessionResult.dirAtoB) == n_windows && isequal(sessionResult.time_centers, time_centers)
                all_dirAtoB(:, sIdx) = sessionResult.dirAtoB;
                all_dirBtoA(:, sIdx) = sessionResult.dirBtoA;
                all_catAtoB(:, sIdx) = sessionResult.catAtoB;
                all_catBtoA(:, sIdx) = sessionResult.catBtoA;
                sessions_processed_count = sessions_processed_count + 1;
                fprintf('  [INFO] Successfully loaded and integrated results for session %d.\n', thisDate);
                continue; % Skip to the next session
            else
                fprintf('  [WARNING] Existing results file %s is incompatible (e.g., different number of windows). Recomputing session.\n', session_mat_filename);
            end
        catch ME_load
            fprintf('  [WARNING] Error loading %s: %s. Recomputing session.\n', session_mat_filename, ME_load.message);
        end
    end
    % --- End Check ---

    % (a) Get Trial Info for this session
    mask = (allSessionIDs == thisDate);
    if ~any(mask)
        fprintf('  [WARNING] No Behavioral (Bhv) data found for session date %d. Skipping.\n', thisDate);
        continue;
    end
    session_bhv = data_master.Bhv(mask);
    if ~isfield(session_bhv, 'Trial_info') || isempty(session_bhv.Trial_info)
         fprintf('  [WARNING] No Trial_info found for session date %d. Skipping.\n', thisDate);
         continue;
    end
    trial_info = session_bhv.Trial_info;

    % Filter trial info: remove NaNs in category, etc.
    if ~isfield(trial_info, 'category') || ~isfield(trial_info, 'direction')
         fprintf('  [WARNING] Trial_info missing category or direction field for session %d. Skipping.\n', thisDate);
         continue;
    end
    catVals = [trial_info.category];
    validMask = ~isnan(catVals);
    trial_info = trial_info(validMask);
    if isempty(trial_info)
        fprintf('  [INFO] No trials with valid category labels for session %d. Skipping.\n', thisDate);
        continue;
    end
    dir_labels = [trial_info.direction]'; % [nTrials x 1]
    cat_labels = [trial_info.category]'; % [nTrials x 1]
    nTrials = length(trial_info);

    fprintf('  Session %d: Found %d trials with valid category labels.\n', thisDate, nTrials);
    if nTrials < min_trials_thresh
         fprintf('  [WARNING] Session %d has only %d trials (<%d threshold). Skipping.\n', thisDate, nTrials, min_trials_thresh);
         continue;
    end

    % (b) Find LIP & FEF neurons recorded on this date
    lipMask = (sess_dates_LIP == thisDate);
    fefMask = (sess_dates_FEF == thisDate);
    neurLIP = allNeurLIP(lipMask);
    neurFEF = allNeurFEF(fefMask);
    nLIP = length(neurLIP);
    nFEF = length(neurFEF);
    if nLIP == 0 || nFEF == 0
        fprintf('  [INFO] No simultaneously recorded LIP (%d units) or FEF (%d units) neurons found for session %d. Skipping.\n', ...
                nLIP, nFEF, thisDate);
        continue;
    end
    fprintf('  Session %d: Found %d LIP neurons and %d FEF neurons.\n', thisDate, nLIP, nFEF);

    % Prepare arrays to store results for THIS session's windows
    this_session_dirAtoB = nan(n_windows, 1);
    this_session_dirBtoA = nan(n_windows, 1);
    this_session_catAtoB = nan(n_windows, 1);
    this_session_catBtoA = nan(n_windows, 1);

    % (c) Loop over time windows for this session
    window_success_flag = false; % Track if any window works
    for w = 1:n_windows
        t_center = time_centers(w);
        t_start  = max(t_center - window_size/2, 0); % Ensure start time >= 0
        t_end    = t_center + window_size/2;
        fprintf('  --- Window %d/%d: center=%.1f ms, [%.1f, %.1f] ms ---\n',...
            w, n_windows, t_center*1e3, t_start*1e3, t_end*1e3);

        % (c.1) Extract spike matrices
        try
            [~, X_LIP_m, Y_FEF_m] = extract_spike_matrices( ...
                trial_info, neurLIP, neurFEF, sample_rate, ...
                [t_start, t_end], alignment_event, sub_bin_size);
        catch ME_extract
             fprintf('    [ERROR extracting spikes @ win %d]: %s. Skipping window.\n', w, ME_extract.message);
             continue;
        end

        % (c.2) Check dimensions and sum over neurons
        min_bins = min(size(X_LIP_m, 2), size(Y_FEF_m, 2));
        if min_bins < 1
            fprintf('    [INFO] Zero time bins extracted for window %d. Skipping window.\n', w);
            continue;
        end
        X_LIP_m = X_LIP_m(:, 1:min_bins, :);
        Y_FEF_m = Y_FEF_m(:, 1:min_bins, :);
        X_pop = sum(X_LIP_m, 3);
        Y_pop = sum(Y_FEF_m, 3);
        [NTr_check, Ntime] = size(X_pop);
        if size(Y_pop, 1) ~= NTr_check || size(Y_pop, 2) ~= Ntime || NTr_check ~= nTrials
             fprintf('    [WARNING @ win %d] Trial dimension mismatch. Skipping window.\n', w);
             continue;
        end
        tauUse = min(FIT_opts_base.tauVal, Ntime - 1);
        if Ntime < 2 || tauUse < 1
            fprintf('    [INFO] Not enough time bins (%d) for tau=%d in window %d. Skipping window.\n', Ntime, FIT_opts_base.tauVal, w);
            continue;
        end

        % (c.3) Reshape data for MINT
        X1 = reshape(X_pop, [1, Ntime, nTrials]);
        Y1 = reshape(Y_pop, [1, Ntime, nTrials]);
        dir_3D = reshape(dir_labels, [1, 1, nTrials]);
        cat_3D = reshape(cat_labels, [1, 1, nTrials]);
         fprintf('    Prepared data shapes for MINT: X1=%s, Y1=%s, dir_3D=%s, cat_3D=%s\n', ...
                 mat2str(size(X1)), mat2str(size(Y1)), mat2str(size(dir_3D)), mat2str(size(cat_3D)));

        % (c.4) Setup FIT options
        FIT_opts = FIT_opts_base;
        FIT_opts.tpres = {Ntime};
        FIT_opts.tau   = {tauUse};
        if FIT_opts.tpres{1} <= FIT_opts.tau{1}
            fprintf('    [WARNING @ win %d] tpres <= tau. Skipping window.\n', w);
            continue;
        end
        outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'};

        % (c.5) Compute FIT for DIRECTION
        fit_dir_success = false;
        try
            inputs_dir = {X1, Y1, dir_3D};
            fprintf('    [Direction] Calling FIT...\n');
            [FIT_vals_dir, ~, ~] = FIT(inputs_dir, outputList, FIT_opts);
            fitAB_dir = mean(FIT_vals_dir{1}, 2);
            fitBA_dir = mean(FIT_vals_dir{2}, 2);
            this_session_dirAtoB(w) = fitAB_dir;
            this_session_dirBtoA(w) = fitBA_dir;
            fprintf('      => DIR-FIT(A->B)=%.5f, DIR-FIT(B->A)=%.5f\n', fitAB_dir, fitBA_dir);
            fit_dir_success = true;
        catch ME_fit_dir
            fprintf('    [ERROR computing direction-FIT @ win %d]: %s\n', w, ME_fit_dir.message);
        end

        % (c.6) Compute FIT for CATEGORY
        fit_cat_success = false;
        try
            inputs_cat = {X1, Y1, cat_3D};
            fprintf('    [Category] Calling FIT...\n');
            [FIT_vals_cat, ~, ~] = FIT(inputs_cat, outputList, FIT_opts);
            fitAB_cat = mean(FIT_vals_cat{1}, 2);
            fitBA_cat = mean(FIT_vals_cat{2}, 2);
            this_session_catAtoB(w) = fitAB_cat;
            this_session_catBtoA(w) = fitBA_cat;
            fprintf('      => CAT-FIT(A->B)=%.5f, CAT-FIT(B->A)=%.5f\n', fitAB_cat, fitBA_cat);
            fit_cat_success = true;
        catch ME_fit_cat
            fprintf('    [ERROR computing category-FIT @ win %d]: %s\n', w, ME_fit_cat.message);
        end

        % Mark session as having some success if either FIT worked
        if fit_dir_success || fit_cat_success
            window_success_flag = true;
        end

    end % End of window loop (w)

    % (d) Process results for THIS session
    if window_success_flag
        % Store in overall matrices
        all_dirAtoB(:, sIdx) = this_session_dirAtoB;
        all_dirBtoA(:, sIdx) = this_session_dirBtoA;
        all_catAtoB(:, sIdx) = this_session_catAtoB;
        all_catBtoA(:, sIdx) = this_session_catBtoA;
        sessions_processed_count = sessions_processed_count + 1;
        fprintf('  --- Finished processing windows for session %d ---\n', thisDate);

        % --- Save individual session results ---
        sessionResult = struct();
        sessionResult.sessionID = thisDate;
        sessionResult.time_centers = time_centers;
        sessionResult.dirAtoB = this_session_dirAtoB;
        sessionResult.dirBtoA = this_session_dirBtoA;
        sessionResult.catAtoB = this_session_catAtoB;
        sessionResult.catBtoA = this_session_catBtoA;
        sessionResult.nTrials = nTrials;
        sessionResult.nLIP = nLIP;
        sessionResult.nFEF = nFEF;
        % Store key parameters used for this session's analysis
        sessionResult.params.alignment_event = alignment_event;
        sessionResult.params.window_size  = window_size;
        sessionResult.params.window_step  = window_step;
        sessionResult.params.sub_bin_size = sub_bin_size;
        sessionResult.params.FIT_opts_base = FIT_opts_base;

        try
            save(session_mat_filename, 'sessionResult');
            fprintf('  [SAVED] Session results to: %s\n', session_mat_filename);
        catch ME_save
            fprintf('  [ERROR] Failed to save session results %s: %s\n', session_mat_filename, ME_save.message);
        end
        % --- End Save individual session results ---

        % --- Plot and Save individual session figures ---
        tAxis = time_centers * 1e3; % Time in ms

        % Direction Plot
        figHandleDir = figure('Name', sprintf('Session %d Direction FIT', thisDate), 'Color', 'w', 'Visible', 'off'); % Keep invisible during save
        hold on;
        plot(tAxis, this_session_dirAtoB, '-r', 'LineWidth', 2, 'DisplayName', 'LIP->FEF (Dir)');
        plot(tAxis, this_session_dirBtoA, '-b', 'LineWidth', 2, 'DisplayName', 'FEF->LIP (Dir)');
        xlabel('Time relative to Stimulus Onset (ms)');
        ylabel('Feature Information Transfer (bits)');
        title(sprintf('Direction FIT - Session %d (N_{tr}=%d, N_{LIP}=%d, N_{FEF}=%d)', thisDate, nTrials, nLIP, nFEF));
        legend('Location', 'best');
        grid on; box off;
        hline = refline(0, 0); hline.Color = 'k'; hline.LineStyle = ':';
        ylim_curr = ylim; ylim([min(-0.001, ylim_curr(1)*1.1), max(0.001, ylim_curr(2)*1.1)]); % Adjust Y-lim slightly

        session_fig_dir_filename = fullfile(output_fig_dir, sprintf('session_%d_FIT_Direction.png', thisDate));
        try
            saveas(figHandleDir, session_fig_dir_filename);
            fprintf('  [SAVED] Session direction plot to: %s\n', session_fig_dir_filename);
        catch ME_fig_save
             fprintf('  [ERROR] Failed to save direction figure %s: %s\n', session_fig_dir_filename, ME_fig_save.message);
        end
        close(figHandleDir); % Close the figure

        % Category Plot
        figHandleCat = figure('Name', sprintf('Session %d Category FIT', thisDate), 'Color', 'w', 'Visible', 'off'); % Keep invisible during save
        hold on;
        plot(tAxis, this_session_catAtoB, '-r', 'LineWidth', 2, 'DisplayName', 'LIP->FEF (Cat)');
        plot(tAxis, this_session_catBtoA, '-b', 'LineWidth', 2, 'DisplayName', 'FEF->LIP (Cat)');
        xlabel('Time relative to Stimulus Onset (ms)');
        ylabel('Feature Information Transfer (bits)');
        title(sprintf('Category FIT - Session %d (N_{tr}=%d, N_{LIP}=%d, N_{FEF}=%d)', thisDate, nTrials, nLIP, nFEF));
        legend('Location', 'best');
        grid on; box off;
        hline = refline(0, 0); hline.Color = 'k'; hline.LineStyle = ':';
        ylim_curr = ylim; ylim([min(-0.001, ylim_curr(1)*1.1), max(0.001, ylim_curr(2)*1.1)]); % Adjust Y-lim slightly

        session_fig_cat_filename = fullfile(output_fig_dir, sprintf('session_%d_FIT_Category.png', thisDate));
         try
            saveas(figHandleCat, session_fig_cat_filename);
            fprintf('  [SAVED] Session category plot to: %s\n', session_fig_cat_filename);
        catch ME_fig_save
             fprintf('  [ERROR] Failed to save category figure %s: %s\n', session_fig_cat_filename, ME_fig_save.message);
        end
        close(figHandleCat); % Close the figure
        % --- End Plot and Save individual session figures ---

    else
        fprintf('  --- No successful windows processed for session %d ---\n', thisDate);
        % Ensure NaN columns remain if session skipped or failed entirely
        all_dirAtoB(:, sIdx) = NaN;
        all_dirBtoA(:, sIdx) = NaN;
        all_catAtoB(:, sIdx) = NaN;
        all_catBtoA(:, sIdx) = NaN;
    end

end % End of session loop (sIdx)

fprintf('\n=== Finished processing/loading all sessions. %d sessions had usable data processed/loaded. ===\n', sessions_processed_count);

%% 6) Calculate and Save Session Averages
if sessions_processed_count > 0
    mean_dirAtoB = nanmean(all_dirAtoB, 2);
    mean_dirBtoA = nanmean(all_dirBtoA, 2);
    mean_catAtoB = nanmean(all_catAtoB, 2);
    mean_catBtoA = nanmean(all_catBtoA, 2);

    avgResults = struct();
    avgResults.time_centers  = time_centers;
    avgResults.dirAtoB_mean = mean_dirAtoB;
    avgResults.dirBtoA_mean = mean_dirBtoA;
    avgResults.catAtoB_mean = mean_catAtoB;
    avgResults.catBtoA_mean = mean_catBtoA;
    avgResults.dirAtoB_all  = all_dirAtoB; % Store all session results too
    avgResults.dirBtoA_all  = all_dirBtoA;
    avgResults.catAtoB_all  = all_catAtoB;
    avgResults.catBtoA_all  = all_catBtoA;
    avgResults.sessions_processed_or_loaded = uniqueSessions(~all(isnan(all_dirAtoB) & isnan(all_dirBtoA) & isnan(all_catAtoB) & isnan(all_catBtoA), 1)); % List sessions that contributed
    avgResults.n_sessions_processed_or_loaded = sessions_processed_count;
    % Store key parameters
    avgResults.params.alignment_event = alignment_event;
    avgResults.params.window_size  = window_size;
    avgResults.params.window_step  = window_step;
    avgResults.params.sub_bin_size = sub_bin_size;
    avgResults.params.FIT_opts_base = FIT_opts_base; % Store base options used

    results_filename = sprintf('FIT_allSessions_avg_corrected_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    try
        save(results_filename, 'avgResults');
        fprintf('\n[SAVED] Averaged results to: %s\n', results_filename);
    catch ME_save
        fprintf('\n[ERROR] Failed to save averaged results: %s\n', ME_save.message);
    end

    %% 7) Plot Averaged Results
    tAxis = time_centers * 1e3; % Time in ms
    N_avg = sessions_processed_count; % Use count of successfully processed/loaded sessions

    % Plot Average Direction FIT
    figAvgDir = figure('Name', 'Average Direction FIT Over Time (Corrected)', 'Color', 'w');
    hold on;
    plot(tAxis, mean_dirAtoB, '-r', 'LineWidth', 2, 'DisplayName', sprintf('LIP->FEF (Dir) (N=%d)', N_avg));
    plot(tAxis, mean_dirBtoA, '-b', 'LineWidth', 2, 'DisplayName', sprintf('FEF->LIP (Dir) (N=%d)', N_avg));
    xlabel('Time relative to Stimulus Onset (ms)');
    ylabel('Average Feature Information Transfer (bits)');
    title(sprintf('Average Direction FIT (N=%d Sessions)', N_avg));
    legend('Location', 'best');
    grid on; box off;
    hline = refline(0, 0); hline.Color = 'k'; hline.LineStyle = ':';
    avg_fig_dir_filename = sprintf('FIT_allSessions_avg_Direction_%s.png', datestr(now,'yyyymmdd_HHMMSS'));
     try
        saveas(figAvgDir, avg_fig_dir_filename);
        fprintf('[SAVED] Average direction plot to: %s\n', avg_fig_dir_filename);
    catch ME_fig_save
         fprintf('[ERROR] Failed to save average direction figure %s: %s\n', avg_fig_dir_filename, ME_fig_save.message);
    end
    % Optional: Keep figure open or close(figAvgDir);


    % Plot Average Category FIT
    figAvgCat = figure('Name', 'Average Category FIT Over Time (Corrected)', 'Color', 'w');
    hold on;
    plot(tAxis, mean_catAtoB, '-r', 'LineWidth', 2, 'DisplayName', sprintf('LIP->FEF (Cat) (N=%d)', N_avg));
    plot(tAxis, mean_catBtoA, '-b', 'LineWidth', 2, 'DisplayName', sprintf('FEF->LIP (Cat) (N=%d)', N_avg));
    xlabel('Time relative to Stimulus Onset (ms)');
    ylabel('Average Feature Information Transfer (bits)');
    title(sprintf('Average Category FIT (N=%d Sessions)', N_avg));
    legend('Location', 'best');
    grid on; box off;
    hline = refline(0, 0); hline.Color = 'k'; hline.LineStyle = ':';
    avg_fig_cat_filename = sprintf('FIT_allSessions_avg_Category_%s.png', datestr(now,'yyyymmdd_HHMMSS'));
     try
        saveas(figAvgCat, avg_fig_cat_filename);
        fprintf('[SAVED] Average category plot to: %s\n', avg_fig_cat_filename);
    catch ME_fig_save
         fprintf('[ERROR] Failed to save average category figure %s: %s\n', avg_fig_cat_filename, ME_fig_save.message);
    end
    % Optional: Keep figure open or close(figAvgCat);

else
    fprintf('\n[INFO] No sessions yielded valid data. No results averaged or plotted.\n');
end

fprintf('\n=== CORRECTED MULTI-SESSION FIT ANALYSIS COMPLETE ===\n');

end

%% Helper function (same as before)
function sess_dates = parseSessionDates(NeurIDs)
    % Parses YYYYMMDD from NeuronIDs like 'YYYYMMDD_ChXX_...'
    sess_dates = nan(size(NeurIDs));
    pat = '^(\d{8})_'; % Looks for 8 digits at the start, followed by _
    for i = 1:numel(NeurIDs)
        if ischar(NeurIDs{i}) || isstring(NeurIDs{i})
            t = regexp(NeurIDs{i}, pat, 'tokens', 'once');
            if ~isempty(t)
                try
                    sess_dates(i) = str2double(t{1});
                catch
                    % Handle cases where conversion fails if needed
                end
            end
        end
    end
end