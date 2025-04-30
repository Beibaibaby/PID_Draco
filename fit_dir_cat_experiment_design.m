function fit_dir_cat_experiment_design()
    % FIT_DIR_CAT_EXPERIMENT_DESIGN
    %
    % Demonstrates how to systematically vary certain FIT parameters (e.g.,
    % window size, binning, tau, etc.) to see how the results
    % for direction-specific and category-specific FIT might change.
    %
    % For each "parameter set," we re-run the code that extracts data and
    % computes:
    %   1) LIP->FEF & FEF->LIP FIT about "direction"
    %   2) LIP->FEF & FEF->LIP FIT about "category"
    %
    % and store them in arrays. Then we plot them all for easy comparison.
    %
    % Requirements:
    %   - data_master loaded (RCT_master_dataset_both_monkeys.mat)
    %   - extract_spike_matrices, preprocess_trial_info on path
    %   - MINT Toolbox with FIT() function on path

    clearvars -except data_master
    close all

    %% Debug Print 1
    fprintf('Starting parameter-experiment script for DIR & CAT FIT (PLUGIN ONLY)...\n');

    %% (1) Load data_master if needed
    if ~exist('data_master','var')
        load('RCT_master_dataset_both_monkeys.mat');
        fprintf('[DEBUG] Loaded data_master from file.\n');
    else
        fprintf('[DEBUG] data_master already in workspace.\n');
    end

    %% (2) Select session & trials
    alignment_event = 'Align_to_cat_stim_on';
    sample_rate     = 40000;   % 40 kHz

    session_ids   = [data_master.Bhv.session_id];
    area_LIP      = 'MLIP';
    area_FEF      = 'MFEF';
    selected_date = 20201211;  % Example session
    trial_info    = data_master.Bhv(selected_date == session_ids).Trial_info;

    % Filter trials: correct-only, aligned
    params             = struct();
    params.alignment   = alignment_event;
    params.correct_only= 1;
    [vert_trials, horizon_trials] = preprocess_trial_info(trial_info, params);

    % Choose which set of trials to analyze. E.g. vertical
    trials = vert_trials;

    % Extract neurons from data_master
    neur_LIP = data_master.Neuro.(area_LIP)( ...
        contains({data_master.Neuro.(area_LIP).NeuronID}, num2str(selected_date)));
    neur_FEF = data_master.Neuro.(area_FEF)( ...
        contains({data_master.Neuro.(area_FEF).NeuronID}, num2str(selected_date)));

    % The labels for direction & category
    dir_labels = [trials.direction]';   % shape [nTrials x 1]
    cat_labels = [trials.category]';    % shape [nTrials x 1] => ±1

    fprintf('[DEBUG] #Trials = %d\n', length(trials));
    fprintf('[DEBUG] LIP neurons = %d, FEF neurons = %d\n', length(neur_LIP), length(neur_FEF));

    %% (3) Define a list of parameter sets to test
    % All sets here use 'bias'='plugin'. We vary window size, sub_bin_size, tau, etc.
    % Use double braces {{ }} for 'bin_method' and 'n_bins' to avoid dimension mismatch
    paramSets = {
        struct('window_size',0.050,'window_step',0.010,'sub_bin_size',0.005, ...
               'bin_method',{{'eqpop','eqpop','none'}}, 'n_bins',{{3,3}}, ...
               'bias','plugin','tau',1, 'label','A: 50ms, sub5ms, tau=1'),
        struct('window_size',0.080,'window_step',0.010,'sub_bin_size',0.005, ...
               'bin_method',{{'eqpop','eqpop','none'}}, 'n_bins',{{4,4}}, ...
               'bias','plugin','tau',1, 'label','B: 80ms, sub5ms, tau=1'),
        struct('window_size',0.050,'window_step',0.010,'sub_bin_size',0.0025, ...
               'bin_method',{{'eqpop','eqpop','none'}}, 'n_bins',{{3,3}}, ...
               'bias','plugin','tau',1, 'label','C: 50ms, sub2.5ms, tau=1'),
        struct('window_size',0.050,'window_step',0.010,'sub_bin_size',0.005, ...
               'bin_method',{{'eqpop','eqpop','none'}}, 'n_bins',{{3,3}}, ...
               'bias','plugin','tau',2, 'label','D: 50ms, sub5ms, tau=2')
    };

    nSets = numel(paramSets);

    % We'll store up to 50 windows worth of data. If a param set uses fewer windows, that's fine.
    maxWindows = 50;
    all_dirAtoB = nan(maxWindows, nSets);
    all_dirBtoA = nan(maxWindows, nSets);
    all_catAtoB = nan(maxWindows, nSets);
    all_catBtoA = nan(maxWindows, nSets);

    time_centers = [];

    %% (4) Loop over each paramSet
    for p = 1:nSets
        pars = paramSets{p};

        fprintf('\n[DEBUG] ParamSet %d/%d: "%s"\n', p, nSets, pars.label);
        fprintf('  window_size=%.3f, step=%.3f, sub_bin=%.4f, bins=(%d,%d), bias=%s, tau=%d\n',...
            pars.window_size, pars.window_step, pars.sub_bin_size, ...
            pars.n_bins{1}, pars.n_bins{2}, pars.bias, pars.tau);

        % define local time centers from 0-300 ms
        start_time = 0.0;
        end_time   = 0.3;
        window_size= pars.window_size; 
        window_step= pars.window_step;
        time_centers_local = start_time : window_step : end_time;
        n_windows_local    = length(time_centers_local);

        if p==1
            time_centers = time_centers_local; % for final plotting
        end

        % Initialize arrays for this paramSet
        FIT_dirAtoB = nan(n_windows_local,1);
        FIT_dirBtoA = nan(n_windows_local,1);
        FIT_catAtoB = nan(n_windows_local,1);
        FIT_catBtoA = nan(n_windows_local,1);

        for idx = 1:n_windows_local
            t_center = time_centers_local(idx);
            t_start  = max(t_center - window_size/2, 0);
            t_end    = t_center + window_size/2;

            % Step (a) Extract spiking => [nTrials x nTimeBins x nNeurons]
            [~, X_LIP, Y_FEF] = extract_spike_matrices( ...
                trials, neur_LIP, neur_FEF, sample_rate, ...
                [t_start, t_end], alignment_event, pars.sub_bin_size);

            % Step (b) same # sub-bins
            min_bins = min(size(X_LIP,2), size(Y_FEF,2));
            if min_bins < 1
                continue
            end
            X_LIP = X_LIP(:,1:min_bins,:);
            Y_FEF = Y_FEF(:,1:min_bins,:);

            % Step (c) sum across neurons
            X_LIP = sum(X_LIP,3); % => [nTrials x nTimeBins]
            Y_FEF = sum(Y_FEF,3);
            [nTrials, nTimeBins] = size(X_LIP);

            if size(Y_FEF,1)~=nTrials || size(Y_FEF,2)~=nTimeBins
                continue
            end
            if nTrials < 30 || nTimeBins < (pars.tau+1)
                continue
            end

            % Step (d) reshape => [1 x nTimeBins x nTrials]
            X1 = reshape(X_LIP, [1,nTimeBins,nTrials]);
            Y1 = reshape(Y_FEF, [1,nTimeBins,nTrials]);

            dir_3D = reshape(dir_labels,[1,1,nTrials]);
            cat_3D = reshape(cat_labels,[1,1,nTrials]);

            % Step (e) Build FIT_opts
            FIT_opts = struct();
            FIT_opts.tpres = {nTimeBins};
            % tau can't exceed nTimeBins-1
            tauVal = min(pars.tau, nTimeBins-1);
            FIT_opts.tau   = {tauVal};

            if FIT_opts.tpres{1} <= FIT_opts.tau{1}
                continue
            end
            FIT_opts.redundancy_measure = 'I_min';
            FIT_opts.bin_method         = pars.bin_method; 
            FIT_opts.n_bins             = pars.n_bins;
            FIT_opts.bias               = pars.bias;   % always 'plugin' here
            FIT_opts.pid_constrained    = true;
            FIT_opts.supressWarnings    = false;
            FIT_opts.computeNulldist    = false;
            FIT_opts.parallel_sampling  = true;
            FIT_opts.xtrp               = 10;
            FIT_opts.shuff              = 20;

            % Direction-FIT
            try
                inputs_dir = {X1, Y1, dir_3D};
                outList    = {'FIT(A->B;S)', 'FIT(B->A;S)'};
                [vals_dir, ~, ~] = FIT(inputs_dir, outList, FIT_opts);

                fitAB_dir = mean(vals_dir{1},2);  % average across time dimension
                fitBA_dir = mean(vals_dir{2},2);
                FIT_dirAtoB(idx) = fitAB_dir;
                FIT_dirBtoA(idx) = fitBA_dir;

            catch ME
                fprintf('[DEBUG] Error in direction-FIT at window %d: %s\n', idx, ME.message);
                continue
            end

            % Category-FIT
            try
                inputs_cat = {X1, Y1, cat_3D};
                [vals_cat, ~, ~] = FIT(inputs_cat, outList, FIT_opts);
                fitAB_cat = mean(vals_cat{1},2);
                fitBA_cat = mean(vals_cat{2},2);
                FIT_catAtoB(idx) = fitAB_cat;
                FIT_catBtoA(idx) = fitBA_cat;

            catch ME
                fprintf('[DEBUG] Error in category-FIT at window %d: %s\n', idx, ME.message);
                continue
            end
        end

        % Store in big arrays
        all_dirAtoB(1:n_windows_local, p) = FIT_dirAtoB;
        all_dirBtoA(1:n_windows_local, p) = FIT_dirBtoA;
        all_catAtoB(1:n_windows_local, p) = FIT_catAtoB;
        all_catBtoA(1:n_windows_local, p) = FIT_catBtoA;
    end

    %% (5) Plot comparison: Direction-based
    figure('Name','Direction-based FIT: Param Comparisons (plugin)','Color','w');
    hold on;
    if isempty(time_centers)
        time_centers = 0:0.01:0.3;  % fallback
    end
    tAxis = time_centers*1e3;
    for p = 1:nSets
        plot(tAxis, all_dirAtoB(1:length(tAxis),p), '-','LineWidth',2, ...
            'DisplayName',['A->B: ', paramSets{p}.label]);
        plot(tAxis, all_dirBtoA(1:length(tAxis),p), '--','LineWidth',2, ...
            'DisplayName',['B->A: ', paramSets{p}.label]);
    end
    xlabel('Time (ms)');
    ylabel('Direction-FIT (bits)');
    legend('Location','best');
    grid on; box off;ye
    title('Direction-FIT Over Time for Different Param Sets (plugin)');
    ylim([-0.005, 0.01]);

    %% (6) Plot comparison: Category-based
    figure('Name','Category-based FIT: Param Comparisons (plugin)','Color','w');
    hold on;
    for p = 1:nSets
        plot(tAxis, all_catAtoB(1:length(tAxis),p), '-','LineWidth',2, ...
            'DisplayName',['A->B: ', paramSets{p}.label]);
        plot(tAxis, all_catBtoA(1:length(tAxis),p), '--','LineWidth',2, ...
            'DisplayName',['B->A: ', paramSets{p}.label]);
    end
    xlabel('Time (ms)');
    ylabel('Category-FIT (bits)');
    legend('Location','best');
    grid on; box off;
    title('Category-FIT Over Time for Different Param Sets (plugin)');
    ylim([-0.005, 0.01]);

    fprintf('\nDone! Plotted direction & category FIT for multiple parameter sets (plugin only).\n');
end