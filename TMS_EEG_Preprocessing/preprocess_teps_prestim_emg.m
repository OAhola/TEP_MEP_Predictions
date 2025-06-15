clear;
close all;
%% start diary and define parameters
dtime = string(datetime);
diary_name = string(strcat('preprocessing_diary_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
addpath("eeglab2024.2\")
eeglab nogui;
%rejection threshold shift for frontal channels and thresholds for rejection based on peaks
frontal_threshold_shift = [1.5 0.75]
pre_peak_threshes_large = [70 1 3]
pre_peak_threshes_small = [2 4 0]
%sound paramaters
sound_lambda = 0.05
sound_n_iters = 5
%highpass filter that is applied before ica and before bad channel
%detection
pre_highpass = 2 %Hz
EEG_post_baseline_range = [-20 -10] %post-stimulus baseline range in ms
post_range = [-0.02 0.3] %range of the post-stimulus data
post_range_true = [0.01 0.15] %more true post-stim range for detecting bad channels and trials
post_range_final = [-0.02 0.15] %final data range of the post-stimulus-focused period
pre_range_init = [-1.25 -0.025] %initial pre-stim range
pre_range_final = [-1.045 -0.045] %final post-stim range
time_window_emg = [-0.1 0.1]
emg_preact_window = [-0.1 -0.01]
ica_prob_threshold = 0.75 %this must be exceeded to reject eye comp
n_ica_components = 35 %number of ica components to compute
ica_plot_window = [5 9] %size of ica component plotting window
rejchan_threshold_pre = 2.5 %z-score threshold for rejecting channels from pre-stim
rejchan_threshold_post = 3 %z-score threshold for rejecting channels from post-stim
sspsir_percentile = 90 %sspsir_percentile
%windows for pulse artifact interpolation
pulse_artifact_window1 = [-5 6]
interpolation_window1 = [1 1]
pulse_artifact_window2 = [-5 8]
interpolation_window2 = [1 1]
filter_order = 4 %filter order given to pop_tesa_filtbutter(), effective filter order of two-sided filter
resample_to = 1000 %sampling freq to resample the pre-stim data to
powerline_range = [48 52] %range to filter out powerline noise (applied only to pre-stim)
filter_range = [2 90] %range to filter the pre-stim data to
%trial rejection params
rejtrial_threshold_global_pre = 1.5
rejtrial_thresholds_local_pre = [6 5]
nums_over_rejtrial_thresholds_local_pre = [0 1]
rejtrial_threshold_global_post = 2.5
rejtrial_thresholds_local_post = [6 5]
nums_over_rejtrial_thresholds_local_post = [0 1]
emg_preact_threshold_muv = 50
threshes_extra_emg = emg_preact_threshold_muv + [10 20 30 40 50]
set(0, 'DefaultFigureVisible', 'off')
sites = {'Tuebingen','Aalto'}; %measurement sites
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\' %where to store the data
for site=sites
    %go through sites and define specific paths
    site_char = char(site);
    directory_name_site = fullfile(datapath_base,strcat('Data_',site_char,"\"));
    files_and_folders = dir(directory_name_site);
    is_subfolder = [files_and_folders.isdir];
    folders = files_and_folders(is_subfolder);
    names = {folders.name};
    subject_names = names(contains(names,"sub"));
    for index = 1:length(subject_names)
        reftep_subject = char(subject_names(index))
        %% define file paths and directories for the subject
        directory_path = char(fullfile(directory_name_site,reftep_subject,"\")); %the directory of the EEG data to be loaded
        filename_tmseeg = char(strcat(reftep_subject,'_task-tep_all_eeg.set')); %.set file name of the eeg data
        %% Load data, and divide EEG and EMG:
        % Divide the data into epochs
        EEG_and_EMG_epoched = pop_loadset(filename_tmseeg,directory_path); %load the data
        channel_names_original = {EEG_and_EMG_epoched.chanlocs.labels}; %all channel names
        channel_names_original_lower = lower(channel_names_original); %all channel names lowered
        %% get EMG data
        emg_channels = channel_names_original(contains(channel_names_original_lower,'emg') | ...
            contains(channel_names_original_lower,'fdi') | contains(channel_names_original_lower,'apb')); %select emg channel names
        if length(emg_channels) ~=2
            disp("more or less than 2 emg channels detected")
        end
        EMG = pop_select(EEG_and_EMG_epoched, 'channel', emg_channels); %select only EMG channels
        for i=1:EMG.nbchan
            EMG.chanlocs(i).type = 'EMG'; %set EMG channel type
        end
        %% get EEG data
        eeg_channels = channel_names_original(~contains(channel_names_original_lower,'emg') & ...
            ~contains(channel_names_original_lower,'fdi') & ~contains(channel_names_original_lower,'apb')); %select eeg channel names
        EEG = pop_select(EEG_and_EMG_epoched, 'channel', eeg_channels); %select eeg channels
        for i=1:EEG.nbchan
            EEG.chanlocs(i).type = 'EEG'; %set EEG channel type
        end
        clear EEG_and_EMG_epoched %not needed anymore, save data

        %% separate pre- and post-stimulation windows
        EEG_pre = pop_select(EEG, 'time', pre_range_init); %pre-stim eeg
        EEG_post = pop_select(EEG, 'time', post_range); %post-stim eeg
        post_raw_name = [reftep_subject,'_EEG_post_with_all_chans_epoched.set']; %save data with all channels
        pop_saveset(EEG_post, 'filename',post_raw_name,'filepath',directory_path);
        pre_raw_name = [reftep_subject,'_EEG_pre_with_all_chans_epoched.set']; %save data with all channels
        pop_saveset(EEG_pre, 'filename',pre_raw_name,'filepath',directory_path);
        EEG_pre_with_all_channels = EEG_pre; %save for interpolating channels later
        clear EEG %not used anymore, save memory

        %% Take care of pulse artifact
        % Remove Pulse artefact
        EEG_post = pop_tesa_removedata(EEG_post, pulse_artifact_window1);
        % Interpolate the missing pulse data
        EEG_post = pop_tesa_interpdata(EEG_post, 'cubic', interpolation_window1);

        %% Detect bad channels from the pre- and post data
        EEG_post_true = pop_rmbase(EEG_post, EEG_post_baseline_range); %remove baseline and then select time
        EEG_post_true = pop_select(EEG_post_true, 'time', post_range_true); %post-stim eeg
        EEG_pre_filtered = pop_tesa_filtbutter(EEG_pre, pre_highpass, [], filter_order, 'highpass');
        EEG_pre_filtered = pop_select(EEG_pre_filtered, 'time', pre_range_final);
        disp("Searching for bad channels from pre-stim...")
        bad_channels_pre = find_bad_channels(EEG_pre_filtered, rejchan_threshold_pre, frontal_threshold_shift, pre_peak_threshes_large, pre_peak_threshes_small);
        disp("Searching for bad channels from post-stim...")
        bad_channels_post = find_bad_channels(EEG_post_true, rejchan_threshold_post, [0 0], [], []);
        clear EEG_post_true EEG_pre_filtered %not used anymore, save memory
        union_bad_elecs = union(bad_channels_pre, bad_channels_post); %combine bad channels detected from pre and post
        chan_names = {EEG_pre.chanlocs.labels}; %existing EEG channel names
        bad_chans = chan_names(union_bad_elecs) %bad channel names

        %% remove bad channels from both pre and post data and set datas to average reference
        EEG_pre = pop_select(EEG_pre, 'nochannel', bad_chans); %remove bad channels from the data
        EEG_post = pop_select(EEG_post, 'nochannel', bad_chans); %remove bad channels from the data
        EEG_pre = pop_reref(EEG_pre, []); %average reference the EEG data before ICA
        EEG_post = pop_reref(EEG_post, []); %average reference the EEG data before ICA

        %% Fit ica on pre and remove from pre and post
        EEG_pre_filtered = pop_tesa_filtbutter(EEG_pre, pre_highpass, [], filter_order, 'highpass'); %highpass filter the data before ICA
        EEG_pre_filtered = pop_select(EEG_pre_filtered, 'time', pre_range_final);
        EEG_pre_filtered = pop_tesa_pcacompress(EEG_pre_filtered, 'compVal', n_ica_components, 'plot','off' );
        EEG_pre_filtered = pop_reref(EEG_pre_filtered, []); %average reference the EEG data before ICA
        EEG_pre_filtered = pop_tesa_fastica(EEG_pre_filtered, 'approach','symm','g','tanh','stabilization','off'); %run ICA
        EEG_pre_filtered = iclabel(EEG_pre_filtered); %run iclabel
        ic_probabilities = EEG_pre_filtered.etc.ic_classification.ICLabel.classifications; %n_channels x classes (or labels)
        possible_labels = EEG_pre_filtered.etc.ic_classification.ICLabel.classes; %possible IC labels
        %find the IC probabilities of being an 'Eye' component
        index_of_label = find(strcmp(possible_labels, 'Eye')); %index of the 'Eye' label
        %get the component indices that exceed the threshold
        component_indices_over_thresh = find(ic_probabilities(:,index_of_label) > ica_prob_threshold)
        ic_probabilities(:,index_of_label)
        %% set the ICA properties from filtered pre-stim to the unfiltered pre- and post-stimulus data and remove the components
        %set weights to EEG_pre
        EEG_pre.icaweights = EEG_pre_filtered.icaweights;
        EEG_pre.icasphere = EEG_pre_filtered.icasphere;
        EEG_pre.icawinv = EEG_pre_filtered.icawinv;
        EEG_pre.icachansind = EEG_pre_filtered.icachansind;
        %set weights to EEG_post
        EEG_post.icaweights = EEG_pre_filtered.icaweights;
        EEG_post.icasphere = EEG_pre_filtered.icasphere;
        EEG_post.icawinv = EEG_pre_filtered.icawinv;
        EEG_post.icachansind = EEG_pre_filtered.icachansind;
        %plot and save ica figure
        pop_topoplot(EEG_pre_filtered, 0, 1:n_ica_components,'Independent components', ica_plot_window, 0, 'electrodes', 'off');
        comps_str = sprintf('%d ',component_indices_over_thresh);
        sgtitle(['Component inds removed: ',comps_str])
        saveas(gcf,[directory_path,reftep_subject,'_ica_topographies.png']);
        close(gcf);
        clear EEG_pre_filtered %this variable is not used anymore, save memory
        %remove the selected components from EEG_pre and EEG_post
        EEG_pre = pop_subcomp(EEG_pre, component_indices_over_thresh, 0);
        EEG_post = pop_subcomp(EEG_post, component_indices_over_thresh, 0);

        %% Reconstruct bad channels in the pre-stim with spherical interpolation
        EEG_pre = pop_interp(EEG_pre, EEG_pre_with_all_channels.chanlocs, 'spherical'); %interpolate bad channels in pre-stim
        %filter and downsample the pre-stimulus data
        EEG_pre = pop_tesa_filtbutter(EEG_pre, filter_range(1), filter_range(2), filter_order, 'bandpass');
        EEG_pre = pop_tesa_filtbutter(EEG_pre, powerline_range(1), powerline_range(2), filter_order, 'bandstop');
        EEG_pre = pop_resample(EEG_pre, resample_to);
        EEG_pre = pop_select(EEG_pre, 'time', pre_range_final); %crop pre-stim to its final time range
        EEG_pre = pop_reref(EEG_pre,[]);

        %% remove bad trials detected from pre-stim from pre-stim, post-stim, and EMG
        bad_trials_pre = find_bad_trials(EEG_pre, rejtrial_threshold_global_pre, rejtrial_thresholds_local_pre, nums_over_rejtrial_thresholds_local_pre)
        if ~isempty(bad_trials_pre)
            EEG_pre = pop_select(EEG_pre,'notrial',bad_trials_pre);
            EEG_post = pop_select(EEG_post,'notrial',bad_trials_pre);
            EMG = pop_select(EMG,'notrial',bad_trials_pre);
        end

        %% Suppress noise and reconstruct bad channels in the post-stim period with SOUND, run SSP--SIR to suppress muscle artifacts
        EEG_post = pop_rmbase(EEG_post, EEG_post_baseline_range);
        %run SOUND
        EEG_post = pop_tesa_sound(EEG_post, 'lambdaValue', sound_lambda, 'iter', sound_n_iters, ...
            'leadfieldInFile' ,[], 'leadfieldChansFile', [], 'replaceChans', char(fullfile(directory_path,post_raw_name)), 'multipleConds', []);
        %run SSP--SIR
        EEG_post = pop_tesa_sspsir(EEG_post, 'artScale','automatic','PC', {'data',[sspsir_percentile]});
        %take care of pulse artifact now on a widened window (it is highly
        %likely that the previous interpolation of the window interpolated
        %it with artifactual data)
        EEG_post = pop_tesa_removedata(EEG_post, pulse_artifact_window2);
        EEG_post = pop_tesa_interpdata(EEG_post, 'cubic', interpolation_window2);
        EEG_post = pop_rmbase(EEG_post, EEG_post_baseline_range);
        EEG_post = pop_reref(EEG_post,[]);

        %% find bad trials from post-stim and remove them from pre-stim, post-stim, and EMG
        EEG_post_in_post_range = pop_select(EEG_post, 'time', post_range_true);
        bad_trials_post = find_bad_trials(EEG_post_in_post_range, rejtrial_threshold_global_post, rejtrial_thresholds_local_post, nums_over_rejtrial_thresholds_local_post)
        if ~isempty(bad_trials_post)
            EEG_pre = pop_select(EEG_pre,'notrial',bad_trials_post);
            EEG_post = pop_select(EEG_post,'notrial',bad_trials_post);
            EMG = pop_select(EMG,'notrial',bad_trials_post);
        end
        EEG_post = pop_select(EEG_post, 'time', post_range_final);

        %% process EMG data
        EMG = pop_select(EMG, 'time', time_window_emg);
        EMG = pop_rmbase(EMG, emg_preact_window*1000 ,[]);
        EMG = pop_tesa_detrend(EMG, 'linear',emg_preact_window*1000); %detrend the preact window
        epsilon = 1/EMG.srate;
        EMG_pre = pop_select(EMG, 'time', emg_preact_window + [0 epsilon]);
        %get bad emg trials and remove them from pre-stim, post-stim, and emg
        bad_emg_trials = find_bad_emg_trials(EMG_pre,emg_preact_threshold_muv, threshes_extra_emg) %bad emg trials
        if ~isempty(bad_emg_trials)
            EEG_pre = pop_select(EEG_pre,'notrial',bad_emg_trials);
            EEG_post = pop_select(EEG_post,'notrial',bad_emg_trials);
            EMG = pop_select(EMG,'notrial',bad_emg_trials);
        end
        %% plot emg results
        time_inds = find(EMG.times > 10 & EMG.times < 50);
        %create a figure of the emg post stim
        figure;
        for ch_ind=1:EMG.nbchan
            subplot(EMG.nbchan,1,ch_ind)
            hold on
            mean_data = mean(squeeze(EMG.data(ch_ind,time_inds,:)),2);
            std_data = std(squeeze(EMG.data(ch_ind,time_inds,:)),0,2);
            fill([EMG.times(time_inds),fliplr(EMG.times(time_inds))],[mean_data + std_data; flipud(mean_data - std_data)],[0.8, 0.8, 0.8],'EdgeColor','none')
            plot([EMG.times(time_inds)], mean_data)
            xlabel("Time (ms)")
            ylabel("Voltage (uV)")
            legend('std','mean')
            title(string(EMG.chanlocs(ch_ind).labels))
        end
        figname_emg = strcat(reftep_subject,' - emg data');
        sgtitle(figname_emg)
        saveas(gcf,[directory_path figname_emg, '.png']);
        close(gcf);
        EEG_post.setname = 'TEPs_processed';
        EEG_pre.setname = 'prestim_eeg_processed';
        EMG.setname = 'EMG_processed';
        figure; pop_timtopo(EEG_post, [-20 145], []);
        saveas(gcf,[directory_path,reftep_subject,'_tep_topos.png']);
        close(gcf);

        %% do consistency checks
        EMG = eeg_checkset(EMG);
        EEG_pre = eeg_checkset(EEG_pre);
        EEG_post = eeg_checkset(EEG_post);
        fields_to_remove = {'event','eventlatency','eventtype','eventurevent'};
        trimmed_epochinfo_EEG_post = rmfield(EEG_post.epoch,fields_to_remove);
        trimmed_epochinfo_EEG_pre = rmfield(EEG_pre.epoch,fields_to_remove);
        if ~isequal(EMG.epoch, EEG_post.epoch) || ~isequal(trimmed_epochinfo_EEG_post, trimmed_epochinfo_EEG_pre)
            disp("structures have differing epoch information")
            return
        end
        %% save the processed files
        pop_saveset(EMG, 'filename',[reftep_subject,'_EMG_processed.set'],'filepath',directory_path);
        pop_saveset(EEG_pre, 'filename',[reftep_subject,'_EEG_pre_processed.set'],'filepath',directory_path);
        pop_saveset(EEG_post, 'filename',[reftep_subject,'_EEG_post_processed.set'],'filepath',directory_path);
    end
end