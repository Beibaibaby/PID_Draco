function [vert_trial_info, horizon_trial_info] = preprocess_trial_info(decoder_trial_info, params)

    % Filter out all non-RCT trials
    filter_other_trials = isnan([decoder_trial_info.category].');
    decoder_trial_info(any(filter_other_trials,2))=[];

    % Valid trials only
    if params.correct_only
        valid_trial_inds = find([decoder_trial_info.trial_error]==0);
        temp_trial_info = decoder_trial_info(valid_trial_inds);
    else
        valid_trial_inds = find([decoder_trial_info.trial_error]==0 | [decoder_trial_info.trial_error]==6);
        temp_trial_info = decoder_trial_info(valid_trial_inds);
    end

    % Filter trials missing the alignment event
    empty_events = cellfun(@isempty,{temp_trial_info.(params.alignment)});
    temp_trial_info = temp_trial_info(~empty_events);

    % Vertical/horizontal separation
    vert_trial_info = temp_trial_info([temp_trial_info.targets_vert]==1);
    vert_trial_info(isnan(extractfield(vert_trial_info, 'Align_to_cat_stim_on')))=[];

    horizon_trial_info = temp_trial_info([temp_trial_info.targets_vert]==0);
    horizon_trial_info(isnan(extractfield(horizon_trial_info, 'Align_to_cat_stim_on')))=[];

end