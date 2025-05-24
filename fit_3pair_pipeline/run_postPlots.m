% run_postPlots.m
addpath(genpath(pwd));

for tauVal = 1:5
    resultsRoot = fullfile(pwd,'results',sprintf('tau%d',tauVal));
    if ~isfolder(resultsRoot), continue; end
    fprintf('[post] tau=%d  ->  %s\n',tauVal,resultsRoot);
    make_pair_averages_smooth(resultsRoot,5,'movmean',99);
end
quit
