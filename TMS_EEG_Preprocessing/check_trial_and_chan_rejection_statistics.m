clear;
close all;
%% define the subject and start diary
dtime = string(datetime);
diary_name = string(strcat('trial_and_chan_rejection_stats_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
addpath("eeglab2024.2\")
eeglab nogui;
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\';
sites = {'Tuebingen','Aalto'};
trial_drops = [];
chan_drops = [];
for site=sites
    %load data for subjects (raw and processed)
    %detect how many channels and trials have been rejected 
    %go through all the sites and define specific paths
    site_char = char(site);
    directory_name_site = fullfile(datapath_base,strcat('Data_',site_char,"\"));
    files_and_folders = dir(directory_name_site);
    is_subfolder = [files_and_folders.isdir];
    folders = files_and_folders(is_subfolder);
    names = {folders.name};
    subject_names = names(contains(names,"sub"));
    for index = 1:length(subject_names)
        reftep_subject = char(subject_names(index));
        %% define file paths and directories for the subject
        directory_path = char(fullfile(directory_name_site,reftep_subject,"\")); %the directory of the EEG data to be loaded
        all_eeg_channels_name_after_dropping = char(strcat(reftep_subject,'_EEG_post_with_all_chans_epoched.set')); %load data
        EEG_epoched = pop_loadset('filename',all_eeg_channels_name_after_dropping,'filepath',directory_path);
        number_of_trials_original = size(EEG_epoched.data,3);
        number_of_channels_original = size(EEG_epoched.data,1);
        %load the pre-processed data
        eeg_file_preprocessed = char(strcat(reftep_subject,'_EEG_post_processed.set'));
        EEG_preprocessed = pop_loadset(eeg_file_preprocessed, directory_path);
        number_of_trials_preprocessed = size(EEG_preprocessed.data,3);
        number_of_reconstructed_channels = sum(EEG_preprocessed.badC);

        %check the percentage of dropping bad trials and channels
        trial_dropping_rate = (number_of_trials_original - number_of_trials_preprocessed) / number_of_trials_original;
        channel_dropping_rate = number_of_reconstructed_channels / number_of_channels_original;

        trial_drops(end+1) = trial_dropping_rate;
        chan_drops(end+1) =  channel_dropping_rate;
    end
end
disp("mean std trial drops")
mean(trial_drops)*100
std(trial_drops)*100
disp("mean std chan drops")
mean(chan_drops)*100
std(chan_drops)*100
     
