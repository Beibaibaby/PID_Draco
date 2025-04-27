%% Load data
clearvars -except data_master
close all

if ~exist('data_master','var')
    load('../Draco_Oliver/RCT_master_dataset_both_monkeys.mat');
end

addpath(genpath('/path/to/MINT/src'));  % <-- Update if needed

%% Settings
alignment_event = 'Align_to_cat_stim_on';
time_window = [0.1, 0.5]; % 100-500 ms
sample_rate = 40000; 

area_LIP = 'MLIP'; 
area_FEF = 'MFEF'; 

session_ids = [data_master.Bhv.session_id];
neur_info_LIP = data_master.Neuro.(area_LIP);
neur_info_FEF = data_master.Neuro.(area_FEF);

%% Initialize pooling containers
all_X_LIP = [];
all_Y_FEF = [];
all_S = [];

uniq_dates = unique(session_ids);

%% Choose number of neurons to fix
num_neurons_LIP = 5;
num_neurons_FEF = 5;

%% Loop across sessions
for u = 1:length(uniq_dates)
    this_date = uniq_dates(u);

    trial_info = data_master.Bhv(this_date == session_ids).Trial_info;
    neur_LIP = neur_info_LIP(contains({neur_info_LIP.NeuronID}, num2str(this_date)));
    neur_FEF = neur_info_FEF(contains({neur_info_FEF.NeuronID}, num2str(this_date)));

    if isempty(neur_LIP) || isempty(neur_FEF)
        disp(['Skipping date ', num2str(this_date), ' due to missing neurons']);
        continue
    end

    % Check if enough neurons
    if length(neur_LIP) < num_neurons_LIP || length(neur_FEF) < num_neurons_FEF
        disp(['Skipping date ', num2str(this_date), ' due to not enough neurons']);
        continue
    end

    % Preselect fixed neurons
    neur_LIP = neur_LIP(1:num_neurons_LIP);
    neur_FEF = neur_FEF(1:num_neurons_FEF);

    % Preprocess trials
    [vert_trials, ~] = preprocess_trial_info(trial_info, struct('alignment', alignment_event, 'correct_only', 1, 'vert', 1));

    if isempty(vert_trials)
        disp(['Skipping date ', num2str(this_date), ' due to no valid trials']);
        continue
    end

    % Extract spike counts
    try
        [S, X_LIP, Y_FEF] = extract_spike_matrices(vert_trials, neur_LIP, neur_FEF, sample_rate, time_window, alignment_event);
    catch
        disp(['Skipping date ', num2str(this_date), ' due to mismatch or error']);
        continue
    end

    % Pool
    all_X_LIP = cat(3, all_X_LIP, X_LIP);
    all_Y_FEF = cat(3, all_Y_FEF, Y_FEF);
    all_S = [all_S; S];
end

%% Check trial numbers
disp(['Total pooled trials: ', num2str(length(all_S))])

%% Prepare MINT inputs
X1 = permute(all_X_LIP, [3 2 1]);
X2 = permute(all_Y_FEF, [3 2 1]);
S = all_S(:)';

inputs = {X1, X2, S};
outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'};

%% MINT options
FIT_opts = struct();
FIT_opts.tau = {10};              
FIT_opts.tpres = {23};            
FIT_opts.redundancy_measure = 'I_min'; 
FIT_opts.bin_method = {'eqpop', 'eqpop', 'none'}; 
FIT_opts.n_bins = {3, 3};         
FIT_opts.bias = 'shuffSub';       
FIT_opts.xtrp = 10; 
FIT_opts.pid_constrained = true;
FIT_opts.supressWarnings = true;
FIT_opts.computeNulldist = false;  
FIT_opts.n_samples = 20;          
FIT_opts.shuffling = {'AB_C'};    
FIT_opts.dim_shuffle = {'Trials'};
FIT_opts.parallel_sampling = true;

%% Run MINT
[FIT_corrected, FIT_plugin, FIT_nullDist] = FIT(inputs, outputList, FIT_opts);

%% Display
disp('Direction Information Transfer Results (LIP -> FEF and FEF -> LIP):')
disp(FIT_corrected)