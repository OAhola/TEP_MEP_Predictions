function bad_channels = find_bad_channels(EEG, threshold, frontal_threshold_shift, peak_thresh, peak_thresh2)
    data = EEG.data;
    channel_locations = EEG.chanlocs;
    data_reshaped = reshape(data, size(data,1), size(data,2)*size(data,3));
    %calculate standard deviation of channels and z-score
    channel_sds = get_channel_sds(data_reshaped);
    z_scores = zscore(channel_sds);
    %size(z_scores)
    bad_channels = search_bad_channels(z_scores, threshold, channel_locations, frontal_threshold_shift,[]);
    if ~isempty(bad_channels)
        n_iterations = 1;
        while true
            mask = true(size(z_scores));
            mask(bad_channels) = false;
            new_z_scores = zscore(channel_sds(mask));
            bad_channels_new = search_bad_channels(new_z_scores, threshold, channel_locations, frontal_threshold_shift, bad_channels);
            if isempty(bad_channels_new) %bad_channels_new, does it bring anything new now?
                break
            else
                bad_channels = union(bad_channels,bad_channels_new);
            end
            n_iterations = n_iterations + 1;
        end
    end
    if ~isempty(bad_channels)
        formatted_array = strjoin({channel_locations(bad_channels).labels},', ');
        fprintf('Rejected channels based on std z-score thresholding: [%s]\n',formatted_array)
        fprintf('Took n_iterations: %d\n',n_iterations)
    end

    %drop channels based on outlying large peaks
    if ~isempty(peak_thresh)
        bad_channels_peaks = find_bad_channels_peaks(data, peak_thresh,channel_locations);
        bad_channels = union(bad_channels, bad_channels_peaks);
        if ~isempty(bad_channels_peaks)
            formatted_array = strjoin({channel_locations(bad_channels_peaks).labels},', ');
            fprintf('Rejected channels based on high prominence peaks: [%s]\n',formatted_array)
        end
    end
    if ~isempty(peak_thresh2) %channels that are constantly full of small peaks
        bad_channels_peaks = find_bad_channels_peaks(data, peak_thresh2,channel_locations);
        bad_channels = union(bad_channels, bad_channels_peaks);
        if ~isempty(bad_channels_peaks)
            formatted_array = strjoin({channel_locations(bad_channels_peaks).labels},', ');
            fprintf('Rejected channels based on low prominence peaks: [%s]\n',formatted_array)
        end
    end
end

function bad_channels_peaks = find_bad_channels_peaks(data, peak_thresh,chanlocs)
    n_channels = size(data,1);
    n_trials = size(data,3);
    channel_peaks = zeros(n_channels,1); %init array
    mean_channel_peaks_trials = zeros(n_channels,1); %init array
    %parameters for peak detection and channel rejection
    min_prominence = peak_thresh(1);
    std_scale = peak_thresh(2);
    equals_or_exceeds = peak_thresh(3);
    for i=1:n_channels
        channel_peaks_trials = zeros(n_trials,1); %init array
        for j = 1:n_trials
            data_channel_abs = abs(data(i,:,j));
            [pks, ~] = findpeaks(data_channel_abs,'MinPeakProminence', min_prominence);
            channel_peaks_trials(j) = length(pks);
        end
        %channel_peaks_trials
        channel_peaks(i) = sum(channel_peaks_trials); %save the number of peaks
        mean_channel_peaks_trials(i) = mean(channel_peaks_trials);
        %chanlocs(i).labels
        %length(pks)
    end
    %standard deviation and mean of the number of peaks across channels
    peaks_std = std(channel_peaks);
    mean_channel_peaks = mean(channel_peaks);
    %find(channel_peaks > (mean_channel_peaks + peaks_std*std_scale))
    bad_channels_peaks = find(channel_peaks > (mean_channel_peaks + peaks_std*std_scale) & (mean_channel_peaks_trials >= equals_or_exceeds)); %what channels exceed the tresholds
    %channel_peaks
end


function searched_bad_channels = search_bad_channels(z_scores, threshold, channel_locations, frontal_threshold_shift, exclude)
    searched_bad_channels = []; %init array for bad channels
    n_channels = length(z_scores) + length(exclude);
    z_score_channel = 1;
    for i=1:n_channels
        if ~ismember(i,exclude)
            threshold_now = get_adjusted_treshold(channel_locations, i, threshold, frontal_threshold_shift);
            %channel_locations(i).labels
            %z_scores(z_score_channel)
            if abs(z_scores(z_score_channel)) > threshold_now %check if the threshold has been exceeded
                %channel_locations(i).labels
                %z_scores(z_score_channel)
                searched_bad_channels(end+1) = i;
            end
            z_score_channel = z_score_channel + 1;
        %else
            %channel_locations(i).labels
        end
    end
end


function channel_sds = get_channel_sds(data)
    n_channels = size(data,1); %number of channels
    channel_sds = zeros(n_channels,1); %init sd array for all channels
    %calculate standard deviation of channels
    for i = 1:n_channels
        channel_sds(i) = std(data(i,:));
    end
end

function channel_sds_trials = get_channel_trial_sds(data)
    n_channels = size(data,1); %number of channels
    n_trials = size(data,3); %number of trials
    channel_sds_trials = zeros(n_channels,n_trials); %init sd array for all channels
    %calculate standard deviation and means of channels
    for i = 1:n_channels
        for j = 1:n_trials
            channel_sds_trials(i,j) = std(data(i,:,j));
        end
    end
end

function threshold_now = get_adjusted_treshold(channel_locations, i, threshold, frontal_threshold_shift)
    %adjust thresholds based on shift values
    fp_af_shift = frontal_threshold_shift(1);
    f_shift = frontal_threshold_shift(2);
    chan_name = lower(channel_locations(i).labels);
    if contains(chan_name,'fp') || contains(chan_name,'af')
        threshold_now = threshold + fp_af_shift; %higher threshold for fp/af channels
    elseif strcmp(chan_name(1),'f')
        if strcmp(chan_name(2),'f') || ~isnan(str2double(chan_name(2))) || strcmp(chan_name(2),'z')
            threshold_now = threshold + f_shift; %higher threshold for f channels
        else
            threshold_now = threshold;
        end
    else
        threshold_now = threshold;
    end
end