function run_fit_combo(srcArea, dstArea, labelType, tauVal)
% RUN_FIT_COMBO( …, tauVal )  – called by each SLURM array task

% ---------- fixed paths -----------------------------------------------
mintRoot   = '/project/bdoiron/draco/MINT';
dataFile   = '/project/bdoiron/dracoxu/Draco_Oliver/RCT_master_dataset_both_monkeys.mat';

% results/tauX/....
saveRoot = fullfile(pwd,'results',sprintf('tau%d',tauVal));
if ~exist(saveRoot,'dir'); mkdir(saveRoot); end

addpath(genpath(mintRoot));

AREA = struct('LIP','MLIP','FEF','MFEF','SC','MSC');
assert(isfield(AREA,srcArea) && isfield(AREA,dstArea));

Nworkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
if isempty(gcp('nocreate')), parpool('local',Nworkers); end

load(dataFile,'data_master');

fit_pair_label_allSessions( ...
    data_master, ...
    AREA.(srcArea), AREA.(dstArea), ...
    labelType, saveRoot, tauVal);          % <-- pass tauVal

fprintf('[DONE] %s→%s %s  tau=%d\n',srcArea,dstArea,labelType,tauVal);
end
