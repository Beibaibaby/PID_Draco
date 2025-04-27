%% Load data
clearvars -except data_master
close all

if ~exist('data_master','var')
    load('../Draco_Oliver/RCT_master_dataset_both_monkeys.mat');
end

% Add MINT toolbox path if needed
addpath(genpath('/path/to/MINT/src'))  % <- <- <- update this path!

%% Settings
alignment_event = 'Align_to_cat_stim_on';
time_window = [0.1, 0.5]; % 100-300 ms after stimulus onset
sample_rate = 40000; % spikes sampled at 40 kHz

areas = fieldnames(data_master.Neuro);

% Select LIP and FEF areas
area_LIP = 'MLIP'; 
area_FEF = 'MFEF'; 

% Select sessions and neurons
session_ids = [data_master.Bhv.session_id];
neur_info_LIP = data_master.Neuro.(area_LIP);
neur_info_FEF = data_master.Neuro.(area_FEF);

%% Choose a session (you can loop over sessions later if needed)
selected_date = 20201211; % example session date

% Find trials and neurons from selected date
trial_info = data_master.Bhv(selected_date == session_ids).Trial_info;
neur_LIP = neur_info_LIP(contains({neur_info_LIP.NeuronID}, num2str(selected_date)));
neur_FEF = neur_info_FEF(contains({neur_info_FEF.NeuronID}, num2str(selected_date)));

% Clean trials
[vert_trials, ~] = preprocess_trial_info(trial_info, struct('alignment', alignment_event, 'correct_only', 1, 'vert', 1));

%% Extract spike counts
[S, X_LIP, Y_FEF] = extract_spike_matrices(vert_trials, neur_LIP, neur_FEF, sample_rate, time_window, alignment_event);

% Reshape to match MINT input style: [1, timepoints, trials]
X1 = permute(X_LIP, [3 2 1]);
X2 = permute(Y_FEF, [3 2 1]);
S = S(:)'; % S must be a row vector of trial labels

num_trials = size(X1,3);
num_timepoints = size(X1,2);

%% Prepare inputs for MINT
inputs = {X1, X2, S};
outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'}; % Compute both directions

FIT_opts = struct();
FIT_opts.tau = {10};              % 5 ms delay
FIT_opts.tpres = {23};            % 10 ms present time
FIT_opts.redundancy_measure = 'I_min'; 
FIT_opts.bin_method = {'eqpop', 'eqpop', 'none'}; 
FIT_opts.n_bins = {3, 3};         % 3 bins only
FIT_opts.bias = 'shuffSub';       
FIT_opts.xtrp = 10; 
FIT_opts.pid_constrained = true;
FIT_opts.supressWarnings = true;
FIT_opts.computeNulldist = false;  % important for significance
FIT_opts.n_samples = 20;          % 20 shuffle samples
FIT_opts.shuffling = {'AB_C'};    
FIT_opts.dim_shuffle = {'Trials'};
FIT_opts.parallel_sampling = true;

%% Run MINT FIT
[FIT_corrected, FIT_plugin, FIT_nullDist] = FIT(inputs, outputList, FIT_opts);

%% Display results
disp('Direction Information Transfer Results (LIP -> FEF and FEF -> LIP):')
disp(FIT_corrected)