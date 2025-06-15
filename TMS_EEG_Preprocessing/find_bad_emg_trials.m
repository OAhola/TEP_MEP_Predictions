function bad_emg_trials = find_bad_emg_trials(EMG_pre,emg_preact_threshold_muv, threshes_extra)
        n_trials = size(EMG_pre.data,3);
        n_channels = size(EMG_pre.data,1);
        bad_emg_trials = [];
        %go through all trials
        for i=1:n_trials
            for chan_number=1:n_channels
                data_now = EMG_pre.data(chan_number,:,i);
                diff = max(data_now) - min(data_now);
                if diff > emg_preact_threshold_muv
                %check if the threshold is exceeded within this trial
                    bad_emg_trials(end+1) = i;
                end
            end
        end
        bad_emg_trials = unique(bad_emg_trials);
        n_bad_trials = length(bad_emg_trials);
        fprintf("Found %d bad emg trials\n",length(bad_emg_trials))
        thresh_ind = 1;
        %check how many trials are rejected
        while thresh_ind <= length(threshes_extra)
            current_threshold = threshes_extra(thresh_ind);
            if n_bad_trials > n_trials/3
                fprintf("Adjusting threshold to %d as %d bad trials were found...\n",current_threshold,n_bad_trials)
                bad_emg_trials = [];
                %go through all trials
                for i=1:n_trials
                    for chan_number=1:n_channels
                        data_now = EMG_pre.data(chan_number,:,i);
                        diff = max(data_now) - min(data_now); %min max diff
                        if diff > current_threshold
                        %check if the threshold is exceeded within this trial
                            bad_emg_trials(end+1) = i;
                        end
                    end
                end
                bad_emg_trials = unique(bad_emg_trials);
            else
                break
            end
            thresh_ind = thresh_ind + 1; %adjust treshold index
        end