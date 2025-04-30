%% pid_over_time_cat_pop.m  –– category-specific FIT (population counts)
clearvars -except data_master
close all

if ~exist('data_master','var')
    load RCT_master_dataset_both_monkeys.mat
end
fprintf('\n=== Category-specific FIT (population counts) ===\n');

%% SETTINGS
alignment_event = 'Align_to_cat_stim_on';
Fs            = 40e3;          % sampling rate (Hz)
t0            = 0.00;          % analysis window start (s)
t1            = 0.30;          % analysis window end   (s)
win           = 0.050;         % 50 ms window
step          = 0.010;         % 10 ms hop
t_centers     = t0:step:t1;
n_win         = numel(t_centers);

area_LIP      = 'MLIP';
area_FEF      = 'MFEF';
session_date  = 20201211;

fit_L2F = nan(n_win,1);
fit_F2L = nan(n_win,1);

%% PULL TRIALS & UNITS
sid  = find([data_master.Bhv.session_id]==session_date,1);
trial_info = data_master.Bhv(sid).Trial_info;

params.alignment   = alignment_event;
params.correct_only = 1;
[vertT,horizT] = preprocess_trial_info(trial_info,params);
trials   = [vertT , horizT];

neur_LIP = data_master.Neuro.(area_LIP)( ...
            contains({data_master.Neuro.(area_LIP).NeuronID}, ...
            num2str(session_date)));
neur_FEF = data_master.Neuro.(area_FEF)( ...
            contains({data_master.Neuro.(area_FEF).NeuronID}, ...
            num2str(session_date)));

fprintf('Session %d   • %d LIP units   • %d FEF units   • %d trials\n', ...
        session_date,numel(neur_LIP),numel(neur_FEF),numel(trials));

%% MAIN LOOP
for w = 1:n_win
    tc   = t_centers(w);
    win_s = [max(tc-win/2,0) , tc+win/2];          % window in seconds
    fprintf('Window %02d/%02d  [%4.0f–%4.0f] ms  ', ...
            w,n_win,win_s(1)*1e3,win_s(2)*1e3);

    % --- spike counts & category labels (per trial)
    [S,Xc,Yc] = extract_spike_counts_label( ...
        trials,neur_LIP,neur_FEF,Fs,win_s,alignment_event,'category');

    if numel(unique(S))<2
        fprintf("(single category) — skipped\n"); continue; end

    nT  = size(Xc,1);                 % # trials
    % collapse units: population spike count
    Lcounts = sum(Xc,2);              % trials×1
    Fcounts = sum(Yc,2);              % trials×1

    % reshape → objects×time×trials  (1×1×nT)
    X1 = reshape(Lcounts.',1,1,nT);
    X2 = reshape(Fcounts.',1,1,nT);
    St = reshape(S          ,1,1,nT);

    if w==1, fprintf('[%d trials]  shapes  X1:%s  X2:%s  S:%s\n', ...
            nT, mat2str(size(X1)),mat2str(s ize(X2)),mat2str(size(St))); end
    if nT<=30
        fprintf("(only %d trials) — skipped\n",nT); continue; end

    % --- FIT options
    opt = struct();
    opt.tpres = {1};          % present bin length
    opt.tau   = {0};          % past offset (0 ⇒ no lag history)
    opt.redundancy_measure = 'I_min';
    opt.bin_method = {'eqpop','eqpop','none'};
    opt.n_bins     = {3,3};
    opt.bias       = 'shuffSub';
    opt.pid_constrained = true;
    opt.supressWarnings = true;
    opt.computeNulldist = false;

    try
        [vals,~,~] = FIT({X1,X2,St}, ...
                         {'FIT(A->B;S)','FIT(B->A;S)'}, opt);
        fit_L2F(w) = vals{1};   fit_F2L(w) = vals{2};
        fprintf("L→F %.4e  F→L %.4e\n",vals{1},vals{2});
    catch ME
        fprintf("⚠️  %s  — skipped\n",ME.message);
    end
end

%% PLOT
figure; hold on;
plot(t_centers*1e3,fit_L2F,'r-','LineWidth',2);
plot(t_centers*1e3,fit_F2L,'b-','LineWidth',2);
xlabel('Time from category-stim onset (ms)');
ylabel('FIT (bits)');
title('Category-specific FIT (population counts)');
legend('LIP → FEF','FEF → LIP','Location','best');
grid on; box off;