function fit_over_time_simple_debug()
    % FIT_OVER_TIME_SIMPLE_DEBUG
    %
    % Computes Feature Information Transfer (FIT) over time between two brain areas
    % (LIP and FEF). It:
    %  1) Summarizes all neurons per area into one channel (to avoid huge dimensional states).
    %  2) Reshapes data so dimension #3 = #trials, dimension #2 = #time (bins).
    %  3) Calls MINT's FIT(A->B;S) and FIT(B->A;S) in each 50 ms window.
    %  4) Averages across sub-bins => one FIT value per window.
    %  5) Adds debug prints to diagnose "different number of trials" issues.
    %
    % Make sure you have:
    %   - data_master loaded
    %   - extract_spike_matrices on path
    %   - preprocess_trial_info on path
    %   - MINT toolbox (with FIT.m) on path

    clearvars -except data_master
    close all

    %% --------------------------------------------------------------------
    % 1) Load data_master if needed
    % ---------------------------------------------------------------------
    if ~exist('data_master','var')
        load('RCT_master_dataset_both_monkeys.mat');
    end
    fprintf('Starting FIT over time analysis (debug mode)...\n');

    %% --------------------------------------------------------------------
    % 2) Basic time-window settings
    % ---------------------------------------------------------------------
    alignment_event = 'Align_to_cat_stim_on';
    sample_rate     = 40000;   % 40 kHz
    start_time      = 0.0;     % seconds
    end_time        = 0.3;     % seconds
    window_size     = 0.050;   % 50 ms window
    window_step     = 0.010;   % 10 ms step
    time_centers    = start_time : window_step : end_time;  % array of window centers
    n_windows       = length(time_centers);

    % Subdivide each 50-ms window into smaller bins for MINT (e.g. 5 ms each)
    sub_bin_size    = 0.005;

    %% --------------------------------------------------------------------
    % 3) Select session & trial info
    % ---------------------------------------------------------------------
    session_ids   = [data_master.Bhv.session_id];
    area_LIP      = 'MLIP';
    area_FEF      = 'MFEF';
    selected_date = 20201211;  % Example session date
    trial_info    = data_master.Bhv(selected_date == session_ids).Trial_info;

    % Filter trial info: correct_only, etc.
    params             = struct();
    params.alignment   = alignment_event;
    params.correct_only= 1;
    [vert_trials, horizon_trials] = preprocess_trial_info(trial_info, params);

    % Extract neurons from data_master
    neur_LIP = data_master.Neuro.(area_LIP)( ...
        contains({data_master.Neuro.(area_LIP).NeuronID}, num2str(selected_date)));
    neur_FEF = data_master.Neuro.(area_FEF)( ...
        contains({data_master.Neuro.(area_FEF).NeuronID}, num2str(selected_date)));

    %% --------------------------------------------------------------------
    % 4) Storage for FIT results
    % ---------------------------------------------------------------------
    FIT_AtoB = nan(n_windows,1);
    FIT_BtoA = nan(n_windows,1);

    %% --------------------------------------------------------------------
    % 5) Main loop over time windows
    % ---------------------------------------------------------------------
    for idx = 1:n_windows

        % define the current 50 ms window
        t_center = time_centers(idx);
        t_start  = max(t_center - window_size/2, 0);
        t_end    = t_center + window_size/2;

        fprintf('\n--- Window %d/%d: center=%.1f ms, [%.1f, %.1f] ms ---\n',...
            idx, n_windows, t_center*1e3, t_start*1e3, t_end*1e3);

        %% (a) Extract spikes => [nTrials x nTimeBins x nNeurons]
        [S, X_LIP, Y_FEF] = extract_spike_matrices( ...
            vert_trials, neur_LIP, neur_FEF, sample_rate, ...
            [t_start, t_end], alignment_event, sub_bin_size);

        fprintf('  After extraction:\n');
        disp(['  size(X_LIP)=', mat2str(size(X_LIP)), ...
              ', size(Y_FEF)=', mat2str(size(Y_FEF)), ...
              ', size(S)=', mat2str(size(S))]);

        %% (b) Force same # sub-bins if needed
        min_bins = min(size(X_LIP,2), size(Y_FEF,2));
        X_LIP = X_LIP(:, 1:min_bins, :);
        Y_FEF = Y_FEF(:, 1:min_bins, :);

        %% (c) Sum across neurons => single channel
        %     => shape becomes [nTrials x nTimeBins]
        X_LIP = sum(X_LIP, 3);
        Y_FEF = sum(Y_FEF, 3);

        [nTrials, nTimeBins] = size(X_LIP);
        fprintf('  Summation => shape(X_LIP)=[%d x %d], shape(Y_FEF)=[%d x %d]\n',...
            nTrials, nTimeBins, size(Y_FEF,1), size(Y_FEF,2));

        % Ensure that Y_FEF also has nTrials x nTimeBins
        if size(Y_FEF,1)~=nTrials || size(Y_FEF,2)~=nTimeBins
            fprintf('  Skipping: mismatch in X_LIP vs Y_FEF dimensions.\n');
            continue
        end
        if length(S)~=nTrials
            fprintf('  Skipping: mismatch in #trials between S=%d and X_LIP=%d.\n',...
                length(S), nTrials);
            continue
        end

        %% (d) Reshape so dimension #3 = #trials (MINT expects last dim = #trials).
        %  => final shape for X1, X2: [1 x nTimeBins x nTrials]
        %  => so dimension #2 is time, dimension #3 is trials.
        %   (In MINT: [nDims x nTime x nTrials])
        X1 = reshape(X_LIP, [1, nTimeBins, nTrials]); 
        Y1 = reshape(Y_FEF, [1, nTimeBins, nTrials]);

        % Similarly, S => [1 x 1 x nTrials] if you want to treat S as single time dimension
        S_3D = reshape(S, [1, 1, nTrials]);

        fprintf('  Final shapes => X1:%s, X2:%s, S:%s\n',...
            mat2str(size(X1)), mat2str(size(Y1)), mat2str(size(S_3D)));

        if nTimeBins<1
            fprintf('  Skipping: no sub-bins.\n');
            continue
        end
        if nTrials<30
            fprintf('  Skipping: not enough trials.\n');
            continue
        end

        %% (e) Setup FIT options
        FIT_opts = struct();
        FIT_opts.tpres              = {nTimeBins}; % "present" time
        FIT_opts.tau                = {min(1,nTimeBins-1)};
        FIT_opts.redundancy_measure = 'I_min';
        FIT_opts.bin_method         = {'eqpop','eqpop','none'}; 
        FIT_opts.n_bins             = {3,3}; 
        FIT_opts.bias               = 'plugin';
        FIT_opts.xtrp               = 10;  
        FIT_opts.shuff              = 20;  
        FIT_opts.pid_constrained    = true;
        FIT_opts.supressWarnings    = false;
        FIT_opts.computeNulldist    = false;
        FIT_opts.parallel_sampling  = true;

        % Note that MINT typically sees dimension #1=neurons/dims, #2=time, #3=trials
        % so eqpop binning is done along the "last dimension" by default.
        % If your version supports specifying 'dim_binning', you can do:
        %   FIT_opts.dim_binning = {'Trials','Trials','none'};
        % but it's not strictly necessary if MINT's default is to bin along dim #3.

        if FIT_opts.tpres{1} <= FIT_opts.tau{1}
            fprintf('  Skipping: tpres <= tau.\n');
            continue
        end

        %% (f) Call FIT
        try
            inputs     = {X1, Y1, S_3D};
            outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'};
            fprintf('  Calling FIT now...\n');
            [FIT_vals, ~, ~] = FIT(inputs, outputList, FIT_opts);

            % MINT returns a 1 x nTime array for each measure if dimension #2 is time
            fitAB = FIT_vals{1}; % shape [1 x nTimeBins]
            fitBA = FIT_vals{2}; % shape [1 x nTimeBins]

            % Average across the time dimension (#2) => single scalar
            FIT_AtoB(idx) = mean(fitAB,2);
            FIT_BtoA(idx) = mean(fitBA,2);

            fprintf('  => FIT(A->B)=%.5f, FIT(B->A)=%.5f\n',...
                FIT_AtoB(idx), FIT_BtoA(idx));

        catch ME
            fprintf('  Error in FIT computation: %s\n', ME.message);
            % optionally, do "keyboard" here for debugging
            continue
        end
    end

    %% --------------------------------------------------------------------
    % 6) Plot results
    % ---------------------------------------------------------------------
    fprintf('\nFinished FIT over time analysis.\n');
    figure('Name','FIT Over Time (Debug)','Color','w');
    hold on;
    tAxis = time_centers*1e3;
    plot(tAxis, FIT_AtoB, '-r','LineWidth',2,'DisplayName','LIP->FEF');
    plot(tAxis, FIT_BtoA, '-b','LineWidth',2,'DisplayName','FEF->LIP');
    xlabel('Time (ms)');
    ylabel('Feature Information Transfer (bits)');
    legend('Location','best');
    grid on; box off;
    title('FIT Over Time: LIP & FEF (Summed), dimension #3=trials');
    ylim([-0.005, 0.005]);
end