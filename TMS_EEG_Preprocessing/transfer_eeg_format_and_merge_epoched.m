%%
clear all;
clc;
dtime = string(datetime);
diary_name = string(strcat('diary_change_data_format_and_epoch_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(strcat('the datetime is ', datestr(now, 'dd/mm/yy-HH:MM')))
addpath("eeglab2024.2") %eeglab in the current folder
eeglab nogui; %start eeglab
electrode_location_path = 'standard_1005.elc';
reftep_eeg_excel = "REFTEP_ppEEG_TueAalto.xlsx";
eeg_table = readtable(reftep_eeg_excel); %the table with eeg information
subids = eeg_table.Subject_id;
subids = unique(subids);
%define minimum and maximum ranges for epochs in seconds
epoch_min = -1.3
epoch_max = 0.6
for sub_ind = 1:length(subids)
    sub_id = char(subids{sub_ind}); %the current subject id
    if ismember(sub_id, {'112','119'})
        continue
    end
    if strcmp(sub_id(1),'0')
        site = 'Tuebingen';
        trigger_names = {'A - Stimulation'};
    elseif strcmp(sub_id(1),'1')
        site = 'Aalto';
        trigger_names = {'A - Stimulation', 'B - Stimulation'};
    else
        return
    end
    site_path = strcat('Data_',site);
    eeg_path = char(fullfile('D:\primeprep\construct_data_struct\Raw_data_eeglab',site_path));
    eeg_path_subject = fullfile(eeg_path,['sub-' sub_id ]); %eeg datapath to put merged data into
    mkdir(eeg_path_subject) %make the dir (and the ones before it if needed)
     %% load subject info
    %go through all the blocks
    indices_subjects_now_eeg = strcmp(eeg_table.Subject_id,sub_id);
    eeg_table_sub = eeg_table(indices_subjects_now_eeg,:);
    blocks_all = eeg_table_sub.Block; %all blocks in the structure
    n_tep_blocks = length(blocks_all(contains(blocks_all,'tep'))); %get the number of tep blocks
    %check that the order of the blocks is preserved when reading the
    %data because tep blocks are merged
    %%
    if ~strcmp(blocks_all{1}, 'rest')
        disp("first block should be rest")
        return
    elseif ~strcmp(blocks_all{2}, 'precoil')
        disp("second block should be precoil (resting state with noise and with coil before tms blocks)")
        return
    elseif ~strcmp(blocks_all{end}, 'postcoil')
       disp("last block should be postcoil (resting state with noise and with coil after tms blocks)")
       return
    end
    for tep_block_number=1:n_tep_blocks
        tep_block_identifier = 'tep_run-0' + string(tep_block_number);
        if ~strcmp(blocks_all{2 + tep_block_number}, tep_block_identifier)
            disp("tep blocks are not in order")
            return
        end
    end
    %%
    %go through all blocks
    tep_block_nro = 1; %init tep block number
    for block_ind = 1:length(blocks_all)
        block_str = char(blocks_all{block_ind})  % Extract the string from the cell array
        indices_blocks_now = strcmp(eeg_table_sub.Block,block_str);
        eeg_table_block = eeg_table_sub(indices_blocks_now,:);
        session_num = eeg_table_block.EEG_session_index;
        subject_in_reftep = ['REFTEP_' sub_id]; %like this in tuebingen
        if strcmp(site,'Aalto')
            subject_in_reftep = ['REFTEP' sub_id]; %change for aalto
        end
        session_path = char(fullfile("D:\REFTEP_ALL\Data_raw",strcat("Raw_NeurOne_",site),subject_in_reftep,eeg_table_block.EEG_session_folder));
        EEG = custom_pop_readneurone_tue_noeloc(session_path, session_num);
        channel_names_original = {EEG.chanlocs.labels};
        channel_names_original_lower = lower(channel_names_original);
        channels = channel_names_original(~contains(channel_names_original_lower,'input')); %select proper channels, this one is bad in tuebingen
        EEG = pop_select(EEG, 'channel', channels);
        EEG = pop_chanedit(EEG,'lookup',electrode_location_path,'nosedir','+Y');
        EEG = pop_editset(EEG, 'setname', ['sub-' sub_id '-' block_str]);
        EEG = eeg_checkset(EEG);
        %% save eeg data to eeglab format (resting states separately to merged epochs)
        if contains(block_str,'tep')
            if ~isempty(EEG.event)
                EEG.event = EEG.event(ismember({EEG.event.type},trigger_names)); %select only trigger events
                epochs_name = strcat('epochs_block_',block_str);
                EEG = pop_epoch(EEG, trigger_names, [epoch_min, epoch_max], 'newname', epochs_name, 'epochinfo', 'yes');
                if EEG.trials > 0 & ~isempty(EEG.event) %check that there is still events and trials after epoching
                    n_trials = EEG.trials; %number of trials in the data
                    if tep_block_nro == 1
                        block_identifiers = repmat(tep_block_nro, 1, n_trials); %save trial-to-block information
                        EEG_merged = EEG;
                    else
                        block_identifiers = [block_identifiers, repmat(tep_block_nro, 1, n_trials)]; %update trial-to-block info
                        EEG_merged = pop_mergeset(EEG_merged,EEG);
                    end
                    tep_block_nro = tep_block_nro + 1;
                end
            end
        else
            filetosave = ['sub-' sub_id '_task-' block_str '_eeg.set'];
            pop_saveset(EEG, 'filename',filetosave,'filepath',eeg_path_subject);
        end
    end
    for trial_ind = 1:length(EEG_merged.epoch) %set block identifiers for epochs
        EEG_merged.epoch(trial_ind).block_identifier = block_identifiers(trial_ind);
    end
    EEG_merged = pop_editset(EEG_merged, 'setname', ['sub-' sub_id '-teps_merged']);
    n_epochs = size(EEG_merged.data,3);
    epoch_inds = 1:n_epochs;
    bad_noises = false;
    if strcmp(sub_id,'105') %no noise masking in some trials
        bad_epochs_noisemasking = 280:285;
        epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
        bad_noises = true;
    end
    if strcmp(sub_id,'120') %no noise masking in some trials
        bad_epochs_noisemasking = 278:290;
        epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
        bad_noises = true;
    end
    if strcmp(sub_id,'110') %no noise masking in some trials
        bad_epochs_noisemasking = 580:590;
        epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
        bad_noises = true;
    end
    if strcmp(sub_id,'118') %no noise masking in this trial
        bad_epochs_noisemasking = 286:287;
        epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
        bad_noises = true;
    end
    if bad_noises
        EEG_merged  = pop_select(EEG_merged, 'trial',epoch_inds);
        n_epochs = size(EEG_merged.data,3);
        epoch_inds = 1:n_epochs;
    end
    uniq_events = unique([EEG_merged.event.urevent]);
    if length(uniq_events) ~= length([EEG_merged.event.urevent]) | length(unique([EEG_merged.event.epoch])) ~= length([EEG_merged.event.epoch])
        [unique_vals, ~, idx] = unique([EEG_merged.event.urevent]);
        counts1 = histcounts(idx,numel(unique_vals));
        duplicates_urevents = unique_vals(counts1>1);
        [unique_vals, ~, idx] = unique([EEG_merged.event.epoch]);
        counts2 = histcounts(idx,numel(unique_vals));
        duplicates_epochs = unique_vals(counts2>1);
        if ~isequal(counts1, counts2)
            disp(counts1)
            disp(counts2)
            return
        end
        formatted_array1 = strjoin(string(duplicates_urevents),', ');
        formatted_array2 = strjoin(string(duplicates_epochs),', ');
        fprintf('dropping duplicate urevents/epochs! [%s] / [%s]\n',formatted_array1, formatted_array2)
        eventurevents = [];
        for t=1:size(EEG_merged.epoch,2)
            eventurevents(end+1) = EEG_merged.epoch(t).eventurevent{1};
        end
        contaminated_trials_binary = ismember(eventurevents,duplicates_urevents);
        fprintf('dropping %d epochs because trials were duplicate urevents\n', sum(contaminated_trials_binary))
        EEG_merged = pop_rejepoch(EEG_merged, contaminated_trials_binary, 0);
        n_epochs = size(EEG_merged.data,3);
        epoch_inds = 1:n_epochs;
    end
    n_epochs = length(epoch_inds);
    %select at maximum 1200 epochs
    if n_epochs > 1200 %restrict to maximum 1200 trials
        n_epochs = 1200;
        disp("more than 1200 epochs so restricting to 1200 epochs")
        epoch_inds = epoch_inds(epoch_inds <= n_epochs);
    end
    EEG_merged  = pop_select(EEG_merged, 'trial', epoch_inds);
    EEG_merged = eeg_checkset(EEG_merged);
    filetosave_merged = ['sub-' sub_id '_task-tep_all_eeg.set'];
    pop_saveset(EEG_merged, 'filename',filetosave_merged,'filepath',eeg_path_subject);
    filetosave_block_identifiers = char(fullfile(eeg_path_subject, ['sub-' sub_id '_block_identifiers.mat']));
    block_identifiers_trials = {EEG_merged.epoch.block_identifier};
    save(filetosave_block_identifiers, 'block_identifiers_trials');
end
diary off