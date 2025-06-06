function bad_trials = find_bad_trials(EEG, threshold_global, thresholds_local, nums_over_tresholds_local)
    data = EEG.data;
    trial_sds =  get_trial_sds(data);
    z_scores = zscore(trial_sds);
    %bad trials are those which exceed the boundary limits
    bad_trials_general = find(abs(z_scores) > threshold_global);
    
    %detect local activity
    channel_trial_sds = get_channel_trials_sds(data);
    z_scores2 = zscore(channel_trial_sds);
    n_trials = size(z_scores2,2);
    bad_trials_chans = [];
    n_thresholds = length(thresholds_local);
    if n_thresholds ~= length(nums_over_tresholds_local) %check
        disp("length of thresholds_local and nums_over_tresholds_local must be identical")
        return
    end
    if n_thresholds ~=0
        for i=1:n_trials
            for thresh_ind = 1:n_thresholds
                threshold_local = thresholds_local(thresh_ind);
                num_over_treshold_local = nums_over_tresholds_local(thresh_ind);
                %how many are over threshold_local?
                num_over_thresh = sum(abs(z_scores2(:,i)) > threshold_local);
                if num_over_thresh > num_over_treshold_local
                    %num_over_thresh
                    %num_over_treshold_local
                    %threshold_local
                    %inds = find(abs(z_scores2(:,i)) > threshold_local);
                    %z_scores2(inds,i)
                    %inds
                    %EEG.chanlocs(inds).labels
                    %i
                    bad_trials_chans(end+1) = i;
                end
            end
        end
    end
    %duplicates are removed with union
    bad_trials = union(bad_trials_general, bad_trials_chans);

end

function trial_sds = get_trial_sds(data)
    n_trials = size(data,3);
    trial_sds = zeros(n_trials,1);
    for i = 1:n_trials
        trial_sds(i) = std(data(:,:,i), 0, 'all'); 
    end
end

function channel_trial_sds = get_channel_trials_sds(data)
    n_trials = size(data,3);
    n_channels = size(data, 1);
    channel_trial_sds = zeros(n_channels, n_trials);
    for i = 1:n_channels
        for j = 1:n_trials
            %sd of each channel at each trial
            channel_trial_sds(i,j) = std(data(i,:,j));
        end
    end
end