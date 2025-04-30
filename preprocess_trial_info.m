function [vert_trial_info, horizon_trial_info] = preprocess_trial_info(decoder_trial_info, params)
% PREPROCESS_TRIAL_INFO  Filters out non-RCT trials, errors, missing alignment events, etc.
%
%   Returns two sets of trials: those with vertical targets and those with horizontal.

    % Filter out all trials that aren't the RCT task
    filter_other_trials = isnan([decoder_trial_info.category].');
    decoder_trial_info(any(filter_other_trials,2))=[];

    % Valid trials only
    if params.correct_only
        valid_inds = find([decoder_trial_info.trial_error] == 0);
        temp_info  = decoder_trial_info(valid_inds);
    else
        valid_inds = find([decoder_trial_info.trial_error] == 0 | ...
                          [decoder_trial_info.trial_error] == 6);
        temp_info  = decoder_trial_info(valid_inds);
    end

    % Filter out trials missing the alignment event
    empty_events = cellfun(@isempty, {temp_info.(params.alignment)});
    temp_info    = temp_info(~empty_events);

    % Vertical/horizontal separation
    vert_trial_info    = temp_info([temp_info.targets_vert] == 1);
    vert_trial_info(isnan(extractfield(vert_trial_info, params.alignment))) = [];

    horizon_trial_info = temp_info([temp_info.targets_vert] == 0);
    horizon_trial_info(isnan(extractfield(horizon_trial_info, params.alignment))) = [];

end