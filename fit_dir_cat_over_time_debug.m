function fit_dir_cat_over_time_debug()
    % FIT_DIR_CAT_OVER_TIME_DEBUG
    % 
    % Computes Feature Information Transfer (FIT) over time between two 
    % brain areas (LIP and FEF), for:
    %   1) The "direction" label (trial_info.direction)
    %   2) The "category" label (trial_info.category, ±1)
    %
    % Similar to fit_over_time_simple_debug, but each time window we do 
    % TWO FIT computations: one with direction as S, one with category as S.

    clearvars -except data_master
    close all

    %% 1) Load data_master if needed
    if ~exist('data_master','var')
        load('RCT_master_dataset_both_monkeys.mat');
    end
    fprintf('Starting DIR & CAT FIT over time analysis (debug mode)...\n');

    %% 2) Basic time-window settings
    alignment_event = 'Align_to_cat_stim_on';
    sample_rate     = 40000;   % 40 kHz
    start_time      = 0.0;     % seconds
    end_time        = 0.3;     % seconds
    window_size     = 0.050;   % 50 ms window
    window_step     = 0.010;   % 10 ms step
    time_centers    = start_time : window_step : end_time;
    n_windows       = length(time_centers);

    % Subdivide each 50-ms window into smaller bins for MINT (e.g. 5 ms each)
    sub_bin_size    = 0.005;

    %% 3) Select session & trial info
    session_ids   = [data_master.Bhv.session_id];
    area_LIP      = 'MLIP';
    area_FEF      = 'MFEF';
    selected_date = 20201211;  % Example session date
    %selected_date = 20201216;  % Example session date
    trial_info    = data_master.Bhv(selected_date == session_ids).Trial_info;

    % Filter trial info: correct_only, etc.
    params             = struct();
    params.alignment   = alignment_event;
    params.correct_only= 1;
    [vert_trials, horizon_trials] = preprocess_trial_info(trial_info, params);

    % ---------------------------------------------------------------------
    % CHOOSE which set of trials you want to analyze:
    %  e.g. all 'vert_trials' or 'horizon_trials' or combine them.
    % For example, let's just do vertical:
    trials = vert_trials;
    % Alternatively, trials = [vert_trials, horizon_trials];
    % up to your preference.
    % ---------------------------------------------------------------------

    % Extract neurons from data_master
    neur_LIP = data_master.Neuro.(area_LIP)( ...
        contains({data_master.Neuro.(area_LIP).NeuronID}, num2str(selected_date)));
    neur_FEF = data_master.Neuro.(area_FEF)( ...
        contains({data_master.Neuro.(area_FEF).NeuronID}, num2str(selected_date)));

    %% 4) We'll store 4 arrays: direction-based & category-based
    FIT_dirAtoB = nan(n_windows,1);
    FIT_dirBtoA = nan(n_windows,1);
    FIT_catAtoB = nan(n_windows,1);
    FIT_catBtoA = nan(n_windows,1);

    % Also extract direction & category from these trials:
    dir_labels = [trials.direction]';   % shape [nTrials x 1]
    cat_labels = [trials.category]';    % shape [nTrials x 1] => ±1

    %% 5) Main loop over time windows
    for idx = 1:n_windows
        t_center = time_centers(idx);
        t_start  = max(t_center - window_size/2, 0);
        t_end    = t_center + window_size/2;

        fprintf('\n--- Window %d/%d: center=%.1f ms, [%.1f, %.1f] ms ---\n',...
            idx, n_windows, t_center*1e3, t_start*1e3, t_end*1e3);

        % (a) Extract spikes => [nTrials x nTimeBins x nNeurons]
        %    using the selected trials
        [S_dummy, X_LIP, Y_FEF] = extract_spike_matrices( ...
            trials, neur_LIP, neur_FEF, sample_rate, ...
            [t_start, t_end], alignment_event, sub_bin_size);

        % The 'S_dummy' from extract_spike_matrices is the direction field 
        % by default in your function. We'll ignore it here since we 
        % explicitly define dir_labels & cat_labels above.

        fprintf('  After extraction:\n');
        disp(['  size(X_LIP)=', mat2str(size(X_LIP)), ...
              ', size(Y_FEF)=', mat2str(size(Y_FEF)), ...
              ', #Trials in "trials" = ', num2str(length(trials))]);

        % (b) Force same # sub-bins if needed
        min_bins = min(size(X_LIP,2), size(Y_FEF,2));
        X_LIP = X_LIP(:, 1:min_bins, :);
        Y_FEF = Y_FEF(:, 1:min_bins, :);

        % (c) Sum across neurons => single channel
        X_LIP = sum(X_LIP, 3);  % shape => [nTrials x nTimeBins]
        Y_FEF = sum(Y_FEF, 3);  % shape => [nTrials x nTimeBins]
        [nTrials, nTimeBins] = size(X_LIP);

        fprintf('  Summation => shape(X_LIP)=[%d x %d], shape(Y_FEF)=[%d x %d]\n',...
            nTrials, nTimeBins, size(Y_FEF,1), size(Y_FEF,2));

        % Check dimension match
        if size(Y_FEF,1)~=nTrials || size(Y_FEF,2)~=nTimeBins
            fprintf('  Skipping: mismatch X_LIP vs Y_FEF.\n');
            continue
        end

        % (d) Reshape so dimension #3 = #trials (MINT expects last dim=trials)
        % => [1 x nTimeBins x nTrials]
        X1 = reshape(X_LIP, [1, nTimeBins, nTrials]);
        Y1 = reshape(Y_FEF, [1, nTimeBins, nTrials]);

        % Also reshape the label vectors:
        %   dir_labels => [nTrials x 1] => [1 x 1 x nTrials]
        %   cat_labels => [nTrials x 1] => [1 x 1 x nTrials]
        dir_3D = reshape(dir_labels, [1, 1, nTrials]);
        cat_3D = reshape(cat_labels, [1, 1, nTrials]);

        fprintf('  Final shapes => X1:%s, X2:%s, dir_3D:%s, cat_3D:%s\n',...
            mat2str(size(X1)), mat2str(size(Y1)), ...
            mat2str(size(dir_3D)), mat2str(size(cat_3D)));

        if nTimeBins<1
            fprintf('  Skipping: no sub-bins.\n');
            continue
        end
        if nTrials<30
            fprintf('  Skipping: not enough trials.\n');
            continue
        end

        %% (e) Setup the FIT options
        FIT_opts = struct();
        FIT_opts.tpres              = {nTimeBins};
        FIT_opts.tau                = {min(1, nTimeBins-1)};
        FIT_opts.redundancy_measure = 'I_min';
        FIT_opts.bin_method         = {'eqpop','eqpop','none'};  % last 'none' => no binning for the label
        FIT_opts.n_bins             = {3,3}; 
        FIT_opts.bias               = 'plugin';    % can switch to 'shuffSub' etc.
        FIT_opts.xtrp               = 10;  
        FIT_opts.shuff              = 20;  
        FIT_opts.pid_constrained    = true;
        FIT_opts.supressWarnings    = false;
        FIT_opts.computeNulldist    = true;
        FIT_opts.parallel_sampling  = true;

        if FIT_opts.tpres{1} <= FIT_opts.tau{1}
            fprintf('  Skipping: tpres <= tau.\n');
            continue
        end

        %% (f) Compute FIT for DIRECTION
        try
            inputs_dir  = {X1, Y1, dir_3D};
            outputList  = {'FIT(A->B;S)', 'FIT(B->A;S)'};
            fprintf('  [Direction] Calling FIT...\n');
            [FIT_vals_dir, ~, ~] = FIT(inputs_dir, outputList, FIT_opts);

            % MINT returns shape [1 x nTimeBins], so average over dimension #2
            fitAB_dir = mean(FIT_vals_dir{1},2);
            fitBA_dir = mean(FIT_vals_dir{2},2);

            FIT_dirAtoB(idx) = fitAB_dir;
            FIT_dirBtoA(idx) = fitBA_dir;
            fprintf('    => DIR-FIT(A->B)=%.5f, DIR-FIT(B->A)=%.5f\n',...
                fitAB_dir, fitBA_dir);

        catch ME
            fprintf('  Error computing direction-FIT: %s\n', ME.message);
            continue
        end

        %% (g) Compute FIT for CATEGORY
        try
            inputs_cat = {X1, Y1, cat_3D};
            % same outputList, same opts
            fprintf('  [Category] Calling FIT...\n');
            [FIT_vals_cat, ~, ~] = FIT(inputs_cat, outputList, FIT_opts);

            fitAB_cat = mean(FIT_vals_cat{1},2);
            fitBA_cat = mean(FIT_vals_cat{2},2);

            FIT_catAtoB(idx) = fitAB_cat;
            FIT_catBtoA(idx) = fitBA_cat;
            fprintf('    => CAT-FIT(A->B)=%.5f, CAT-FIT(B->A)=%.5f\n',...
                fitAB_cat, fitBA_cat);

        catch ME
            fprintf('  Error computing category-FIT: %s\n', ME.message);
            continue
        end
    end

    %% 6) Plot results for direction
    fprintf('\nFinished Direction & Category FIT over time analysis.\n');

    figure('Name','Direction-based FIT Over Time','Color','w');
    hold on;
    tAxis = time_centers*1e3;
    plot(tAxis, FIT_dirAtoB, '-r','LineWidth',2,'DisplayName','LIP->FEF (Dir)');
    plot(tAxis, FIT_dirBtoA, '-b','LineWidth',2,'DisplayName','FEF->LIP (Dir)');
    xlabel('Time (ms)');
    ylabel('Feature Information Transfer (bits)');
    legend('Location','best');
    grid on; box off;
    title('FIT Over Time: Direction Labels');
    ylim([-0.005, 0.005]);

    %% 7) Plot results for category
    figure('Name','Category-based FIT Over Time','Color','w');
    hold on;
    plot(tAxis, FIT_catAtoB, '-r','LineWidth',2,'DisplayName','LIP->FEF (Cat)');
    plot(tAxis, FIT_catBtoA, '-b','LineWidth',2,'DisplayName','FEF->LIP (Cat)');
    xlabel('Time (ms)');
    ylabel('Feature Information Transfer (bits)');
    legend('Location','best');
    grid on; box off;
    title('FIT Over Time: Category Labels (±1)');
    ylim([-0.005, 0.005]);

end