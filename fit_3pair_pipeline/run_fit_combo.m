function run_fit_combo(srcArea, dstArea, labelType)
% RUN_FIT_COMBO  – called by each SLURM array task
%
% Examples:
%   run_fit_combo('LIP','FEF','direction')
%   run_fit_combo('FEF','SC','category')

% ---------- fixed paths (edit if you move things) -----------------------
mintRoot   = '/project/bdoiron/draco/MINT';                 % compiled toolbox
dataFile   = '/project/bdoiron/dracoxu/Draco_Oliver/RCT_master_dataset_both_monkeys.mat';
saveRoot   = fullfile(pwd,'results');                       % under pipeline folder
if ~exist(saveRoot,'dir'); mkdir(saveRoot); end

addpath(genpath(mintRoot));

% ---------- area names in data_master -----------------------------------
AREA = struct('LIP','MLIP','FEF','MFEF','SC','MSC');
assert(isfield(AREA,srcArea) && isfield(AREA,dstArea), 'Area names bad');

% ---------- open local parallel pool (match SLURM cpus) -----------------
Nworkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
if isempty(gcp('nocreate'))
    parpool('local',Nworkers);
end

% ---------- load data once ----------------------------------------------
load(dataFile,'data_master');

% ---------- compute ------------------------------------------------------
fit_pair_label_allSessions( ...
    data_master, ...
    AREA.(srcArea), AREA.(dstArea), ...
    labelType, saveRoot);

disp('[DONE] combo finished.');
end
