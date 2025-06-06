clear;
clc;
addpath("eeglab2024.2")
addpath("neurone_tools_for_matlab_1.1.3.11")
eeglab nogui;
dtime = string(datetime);
diary_name = string(strcat('aligning_triggers_diary_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
epoch_min = -1.3
epoch_max = 0.6
epoch_window = [-1.3 0.6]
%read excel
dataframe = readtable("REFTEP_ppEEG_TueAalto.xlsx");
subject_names = dataframe.Subject_id;
subject_names = unique(subject_names);
for sub_ind = 1:length(subject_names)
    subject_name = char(subject_names{sub_ind});
    %determine site-specific parameters (take into account the ITI)
    if ~ismember(subject_name,{'119','112'})
        if strcmp(subject_name(1),"0")
            min_difference_threshold = 1200; %closeness of stimuli to be detected in ms
            site = 'Tuebingen';
            system= 'localite';
            offset_detection_thresh = 1200;
            min_lat_thresh = 2500;
            minutes_threshold = 2; %times used to detect breaks to determine offsets between neurone and neuronavigation systems
            selected_event_types_eeglab = {'A - Stimulation'};
            trigger_info_path = strcat('D:\REFTEP_ALL\Neuronavigation_nexstim_localite\',site,'_localite\stimulation_times\sub-',subject_name,'\sub-',subject_name,'_stimulations.mat');
        elseif strcmp(subject_name(1),"1")
            minutes_threshold = 2; %times used to detect breaks to determine offsets between neurone and neuronavigation systems
            min_difference_threshold = 1800; %closeness of stimuli to be detected
            offset_detection_thresh = 1800;
            min_lat_thresh = 3990; %jitters ok
            site = 'Aalto';
            system = 'nexstim';
            selected_event_types_eeglab = {'A - Stimulation','B - Stimulation'};
            trigger_info_path = strcat('D:\REFTEP_ALL\Neuronavigation_nexstim_localite\',site,'_nexstim\stimulation_times\REFTEP',subject_name,'\REFTEP',subject_name,'_stimulations.mat');
        else
            return
        end
    else
        continue
    end

    %find the stimulation file
    trigger_markers_remaining_path = strcat('D:\REFTEP_ALL\Features\Features_',site,'\sub-',subject_name,'\sub-',subject_name,'_stimulations_final.mat');
    neurone_stimulation_times_remaining_path = strcat('D:\REFTEP_ALL\Features\Features_',site,'\sub-',subject_name,'\sub-',subject_name,'_neurone_stimulation_times.mat');
    selected_event_type_module = 'Stimulation';
    subject_df = dataframe(strcmp(subject_name, subject_names),:);
    try
        trigger_info = load(trigger_info_path);
        found_triggers = true;
    catch
        disp("no trigger markers found")
        found_triggers = false;
    end
    %load preprocessed data with the correct number of urevents left and match them
    processed_eeg_filepath = strcat('D:\REFTEP_ALL\EEG_preprocessing_data\',strcat('Data_',site,'\sub-',subject_name),'\sub-',subject_name, '_EEG_post_processed.set');
    eegs = {};
    %% check if the processed_eeg exists
    try
        EEG_preprocessed = pop_loadset(processed_eeg_filepath);
        found_eeg = true;
    catch
        disp("no processed eeg found")
        found_eeg = false;
    end
    if found_eeg %both eeg and trigger file has been found, load the eeg data in two formats
        session_numbers = [];
        included_epochs = {};
        indices_subjects_now_eeg = strcmp(dataframe.Subject_id,subject_name);
        eeg_table_sub = dataframe(indices_subjects_now_eeg,:);
        blocks_all = eeg_table_sub.Block;
        blocks = blocks_all(contains(blocks_all,'tep')); %tep blocks
        good_set_ind = 1;
        bad_block_inds = [];
        for block_ind = 1:length(blocks)
            block_str = char(blocks{block_ind});  % Extract the string from the cell array
            indices_blocks_now = strcmp(eeg_table_sub.Block,block_str);
            eeg_table_block = eeg_table_sub(indices_blocks_now,:);
            session_num = eeg_table_block.EEG_session_index;
            session_numbers(end+1)=session_num;
            if strcmp(site,'Tuebingen')
                subject_in_reftep = append(['REFTEP_' subject_name]); %tuebingen
            elseif strcmp(site,'Aalto')
                subject_in_reftep = append(['REFTEP' subject_name]); %aalto
            else
                return
            end
            session_path = char(fullfile('D:\REFTEP_ALL\Data_raw',strcat("Raw_NeurOne_",site),subject_in_reftep,eeg_table_block.EEG_session_folder));
            EEG_eeglab = custom_pop_readneurone_tue_noeloc(session_path ,session_num);
            eegs{block_ind} = module_read_neurone(session_path ,session_num);
            [EEG_epoched, eventinds] = pop_epoch(EEG_eeglab, selected_event_types_eeglab, epoch_window, 'newname', 'epochs', 'epochinfo', 'yes');
            [EEG_epoched_small, eventinds2] = pop_epoch(EEG_eeglab, selected_event_types_eeglab, [-0.001, 0.001], 'newname', 'epochs', 'epochinfo', 'yes');
            if isempty(EEG_eeglab.event)
                disp("no events found in EEG_eeglab")
                EEG_epoched
                EEG_epoched_small
                session_num
                eeg_table_block
                event_ids = [];
            else
                event_ids = find(ismember({EEG_eeglab.event.type}, selected_event_types_eeglab));
            end
            %event_ids = find(ismember({EEG_eeglab.event.type}, selected_event_types_eeglab));
            if length(event_ids) ~= length(eventinds2)
                disp("even small epoch range cropping did not work here")
                return
            end
            included_epochs{end+1} = ismember(eventinds2,eventinds); %check which epochs have been dropped
          
            if block_ind == 1
               EEG_merged_full = EEG_eeglab;
            else
               EEG_merged_full = pop_mergeset(EEG_merged_full, EEG_eeglab);
            end
            if EEG_epoched.trials > 0 & ~isempty(EEG_epoched.event) 
                EEG_epoched.event = EEG_epoched.event(ismember({EEG_epoched.event.type},selected_event_types_eeglab)); %select only trigger events
                if good_set_ind == 1
                    EEG_merged = EEG_epoched;
                else
                    EEG_merged = pop_mergeset(EEG_merged,EEG_epoched);
                end
                good_set_ind = good_set_ind + 1;
            end
        end
        %post processing to EEG merged as in the file reader
        EEG_merged = pop_editset(EEG_merged, 'setname', ['sub-' subject_name '-teps_merged']);
        n_epochs = size(EEG_merged.data,3);
        n_epochs_orig = n_epochs;
        epoch_inds = 1:n_epochs;
        bad_noises = false;
        if strcmp(subject_name,'105') %no noise masking in some trials
            bad_epochs_noisemasking = 280:285;
            epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
            bad_noises = true;
        end
        if strcmp(subject_name,'120') %no noise masking in some trials
            bad_epochs_noisemasking = 278:290;
            epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
            bad_noises = true;
        end
        if strcmp(subject_name,'110') %no noise masking in some trials
            bad_epochs_noisemasking = 580:590;
            epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
            bad_noises = true;
        end
        if strcmp(subject_name,'118') %no noise masking in this trial
            bad_epochs_noisemasking = 286:287;
            epoch_inds = get_epoch_inds_good(EEG_merged,bad_epochs_noisemasking);
            bad_noises = true;
        end
        if bad_noises
            EEG_merged  = pop_select(EEG_merged, 'trial',epoch_inds);
            n_epochs = size(EEG_merged.data,3);
            epoch_inds_good = epoch_inds;
            epoch_inds = 1:n_epochs;
        else
            epoch_inds_good = epoch_inds;
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
            epoch_inds_good = epoch_inds(~contaminated_trials_binary);
            epoch_inds = 1:n_epochs;
        end
        n_epochs = length(epoch_inds);
        %select at maximum 1200 epochs
        if n_epochs > 1200 %restrict to maximum 1200 trials
            n_epochs = 1200;
            disp("more than 1200 epochs so restricting to 1200 epochs")
            epoch_inds = epoch_inds(epoch_inds <= n_epochs);
            epoch_inds_good = epoch_inds_good(1:n_epochs);
        end
        EEG_merged  = pop_select(EEG_merged, 'trial', epoch_inds);
        EEG_merged = eeg_checkset(EEG_merged);
        epoched_eeg_filepath = strcat('D:\REFTEP_ALL\EEG_preprocessing_data\',strcat('Data_',site,'\sub-',subject_name),'\sub-',subject_name,'_task-tep_all_eeg.set');
        EEG_epoched_used = pop_loadset(epoched_eeg_filepath);
        if ~isequal(EEG_epoched_used.event,EEG_merged.event)
            disp("not equal events in EEG_merged and EEG_epoched_used")
            return
        end
        %% get the times of the triggers in real world time
        trigger_times = [];
        urevents = [];
        first_stim_time = NaN;
        for i=1:length(eegs)
            EEG_module = eegs{i};
            % get the starting time of the EEG session
            starting_times = {EEG_module.Session.TableSessionPhase.StartDateTime}; %starting times of session phases
            starting_time = starting_times{session_numbers(i)}; %the starting time of this session phase
            %get the starting time in milliseconds in day time
            starting_splitted = split(starting_time, "T");
            hours_minutes_seconds_time_diff = starting_splitted(2);
            hours_minutes_seconds_time_diff_splitted = split(hours_minutes_seconds_time_diff,"+");
            utc_diff = hours_minutes_seconds_time_diff_splitted(2);
            hours_minutes_seconds = hours_minutes_seconds_time_diff_splitted(1);
            starting_hours = str2double(hours_minutes_seconds{1}(1:2));
            starting_minutes = str2double(hours_minutes_seconds{1}(4:5));
            starting_seconds = str2double(hours_minutes_seconds{1}(7:9));
            starting_milliseconds = str2double(hours_minutes_seconds{1}(10));
            starting_time_in_milliseconds = starting_hours*60*60*1000 + starting_minutes*60*1000 + starting_seconds*1000 + starting_milliseconds;
            
            %get the stimulation times
            EEG_stimulation_times = EEG_module.markers.time;
            %select triggers corresponding to the event of interest
            trigger_indices_of_interest = find(strcmp(EEG_module.markers.type,selected_event_type_module));
            trigger_times_of_interest = EEG_module.markers.time(trigger_indices_of_interest);
            trigger_times_of_interest_in_ms = trigger_times_of_interest*1000;
            trigger_times_of_interest_in_real_world = [starting_time_in_milliseconds + trigger_times_of_interest_in_ms];
            %mean_latency = mean(diff(trigger_times_of_interest_in_real_world))
            trigger_times_of_interest_in_real_world_restricted = trigger_times_of_interest_in_real_world(included_epochs{i});
            if isnan(first_stim_time) && ~isempty(trigger_times_of_interest_in_real_world) %save the first stim time
                [first_stim_time, ~] = min(trigger_times_of_interest_in_real_world);
            end
            if length(trigger_times_of_interest_in_real_world) > 1
                mean_latency = mean(diff(trigger_times_of_interest_in_real_world))
                min_lat = min(diff(trigger_times_of_interest_in_real_world))
                %max_lat = max(diff(trigger_times_of_interest_in_real_world))
                std_latency = std(diff(trigger_times_of_interest_in_real_world))
                ntrigs = length(trigger_times_of_interest_in_real_world)
                if min_lat < min_lat_thresh
                    fprintf("failed latency block min diff check %d %d %d\n",i)
                    diff(trigger_times_of_interest_in_real_world)
                end
            end
            if isempty(trigger_times_of_interest_in_real_world)
                fprintf("empty trigger times in block %d",i)
            end
            trigger_times = [trigger_times; trigger_times_of_interest_in_real_world_restricted];
        end
        n_epochs_orig_trigtimes = length(trigger_times);
        trigger_times = trigger_times(epoch_inds_good); %select the same epoch inds
        if n_epochs_orig_trigtimes ~= n_epochs_orig || length(trigger_times) ~= EEG_merged.trials
            disp("mismatcing number of epochs!")
            return
        end
        %% match events
        %get the matching events (urevents) that we are interested in
        urevents_of_interest = [EEG_merged.event.urevent];
        urevents_of_interest_unique = unique(urevents_of_interest);
        if length(urevents_of_interest_unique) ~= length(urevents_of_interest)
            fprintf("NON-unique urevents in EEG_merged: length(urevents_of_interest_unique) ~= length(urevents_of_interest) %d is not %d\n",length(urevents_of_interest_unique),length(urevents_of_interest))
        end

        if length(urevents_of_interest_unique) == length(trigger_times)
            %then create an array with time stamps and trigger times
            urevents_times = struct('times',trigger_times,'urevents',urevents_of_interest_unique, 'timezone',utc_diff);
            %save the struct, commented out so far
            %save(urevents_times_path,'urevents_times')
        else
            disp("Failed: length(urevents_of_interest) ~= length(trigger_times)")
            return
        end
        urevents_remaining = [EEG_preprocessed.event.urevent]; %remaining events in EEG_preprocessed
        potential_urevents = [urevents_times.urevents]; %all potential urevents that were used originally
        if length(urevents_remaining) ~= length(unique(urevents_remaining))
            fprintf("NON-unique urevents in EEG_preprocessed: length(urevents_of_interest_unique) ~= length(urevents_of_interest), %d vs %d\n",length(urevents_remaining),length(unique(urevents_remaining)))
        end
        [is_member3, ~] = ismember(potential_urevents,urevents_remaining);
        ureventstimestimes = urevents_times.times;
        ureventstimesurevents = urevents_times.urevents;
        ureventstimetimezone= urevents_times.timezone;
        urevents_remaining_info = struct();
        urevents_remaining_info.times =  double(ureventstimestimes(is_member3));
        urevents_remaining_info.urevents = ureventstimesurevents(is_member3);
        trial_arr = 1:size(EEG_merged.epoch,2);
        [is_member_epoch, ~] = ismember(potential_urevents,[EEG_preprocessed.epoch.eventurevent]);
        if ~isequal(is_member_epoch,is_member3)
            disp("member arrays not equal!")
            return
        end
        used_trials = trial_arr(is_member3==1);
        urevents_remaining_info.epochs_inds = used_trials;
        urevents_remaining_info.timezone = ureventstimetimezone;
        urevents_remaining_info.first_stim_time = first_stim_time; %time of the first ever stimulation, just in case
        if EEG_preprocessed.trials ~= length(urevents_remaining_info.urevents)
            disp("wrong number of trials compared to events")
            return
        else
            disp("saving stimulation times for all trials")
            %save the stimulation latencies times of the pre-processed data
            if ~all(diff(urevents_remaining_info.times)>0) || ~all(diff(urevents_remaining_info.urevents)>0) || ~all(diff(urevents_remaining_info.epochs_inds)>0)
                disp("not sorted urevents_remaining_info")
                return
            end
            save(neurone_stimulation_times_remaining_path,'urevents_remaining_info')
        end
        if ~found_triggers
            disp("triggers not found, skipping...")
            continue
        end
        %% search for the breaks in neurone and neuronavigation/stimulus system, save the break lengths and the start of new_sequence indices.
        if strcmp(system, 'nexstim') | strcmp(system, 'localite')
            break_lengths_neuronavi = [];
            break_starts_inds_neuronavi = [];
            for i=2:length(trigger_info.trigger_markers)

                difference = double(trigger_info.trigger_markers{i}.recording_time_in_real_world) - double(trigger_info.trigger_markers{i-1}.recording_time_in_real_world);
                if difference >= minutes_threshold*60*1000  %in ms, if this is true then there is a break in stimulation
                    break_lengths_neuronavi(end+1) = difference;
                    break_starts_inds_neuronavi(end+1) = i;
                end
            end
           if isempty(break_starts_inds_neuronavi)
               %there are no break indices so match by nearest time points,
               %this should be reliable with localite at least
               seq_start = double(trigger_info.trigger_markers{1}.recording_time_in_real_world);
               disp("creating offset with the smallest distance between times")
               [offset_found, index_offset] = min(abs(ureventstimestimes - seq_start));
               pot_offsets = seq_start - ureventstimestimes;
               fprintf("potentially bad offset %d",offset_found)
               found_odd_offset = true;
               offsets = [pot_offsets(index_offset)];
           else
               found_odd_offset = false;
               break_lengths_neurone = [];
               break_starts_inds_neurone = [];
               for i=2:length(ureventstimestimes)
                   difference = double(ureventstimestimes(i)) - double(ureventstimestimes(i-1));
                   if difference >= minutes_threshold*60*1000 %in ms, if this is true then there is a break in stimulation
                       break_lengths_neurone(end+1) = difference;
                       break_starts_inds_neurone(end+1) = i;
                   end
               end
               %now the indices of the starts of new blocks have been found,
               %however, there can be some additional ones especially in the neuronavigation logs
               %e.g. due to breaks between mapping and stimulation starts etc.
            
               %choose what breaks to take from the neuronavigation sequence
               break_inds_neuronavi_from_neurone = [];
               min_differences = [];
               bad_inds = [];
               for i=1:length(break_lengths_neurone)
                   break_len = break_lengths_neurone(i);
                   %differences in block lengths
                   differences_break = abs(break_lengths_neuronavi - break_len);
                   [min_diff, argmin] = min(differences_break); %match the blocks
                   neuronavi_index_start = break_starts_inds_neuronavi(argmin);
                   if min_diff >= min_difference_threshold %can not exceed the threshold
                       bad_inds(end+1) = i;
                   end
                   %check if the index has already been used
                   %if the sequence_started has changed then the offset can change
                   if sum(ismember(break_inds_neuronavi_from_neurone,neuronavi_index_start)) > 0 && sum(ismember(bad_inds,i))==0
                       index = find(break_inds_neuronavi_from_neurone==neuronavi_index_start);
                       fprintf("Matched at least one break wrongly with break len diff %f vs %f. Checking differences in starting times.\n",min_diff,min_differences(index))
                       if min_diff < min_differences(index)
                           bad_inds(end+1) = index;
                       else
                           bad_inds(end+1) = i;
                       end
                   end
                   break_inds_neuronavi_from_neurone(end+1) = neuronavi_index_start;
                   min_differences(end+1) = min_diff;
               end
               %define the time offset of the neuronavigation system and the neurone system
               differences = [];
               break_inds_neuronavi_proper = [];
               break_inds_neuronavi = [];
               offsets = [];
               for i=1:length(break_inds_neuronavi_from_neurone)
                   if sum(ismember(bad_inds,i)) == 0
                       if i==1
                           potential_break_ind = break_inds_neuronavi_from_neurone(i)-1;
                           potential_difference = double(trigger_info.trigger_markers{break_inds_neuronavi_from_neurone(i)-1}.recording_time_in_real_world) - double(ureventstimestimes(break_starts_inds_neurone(i)-1));
                           if abs(potential_difference) < min_difference_threshold
                               break_inds_neuronavi_proper(end+1) = potential_break_ind;
                               differences(end+1) = potential_difference;
                               break_inds_neuronavi(end+1) = potential_break_ind;
                               offsets(end+1) = potential_difference;
                           end
                       end
                       break_inds_neuronavi_proper(end+1) = break_inds_neuronavi_from_neurone(i);
                       differences(end+1) = double(trigger_info.trigger_markers{break_inds_neuronavi_from_neurone(i)}.recording_time_in_real_world) - double(ureventstimestimes(break_starts_inds_neurone(i)));
                   end
               end
               if ~all(diff(break_inds_neuronavi_proper)>0)
                   disp("not increasing break inds")
                   return
               end
               fprintf('differences found original = %g\n',differences/1000)
               %check that are there any suspiciously bad differences
               while true
                   diffs = abs(differences' - differences);
                   [bad_rows, bad_cols] = find(diffs>offset_detection_thresh);
                   if length(unique(bad_rows))~=length(bad_rows) && length(unique(bad_cols))~=length(bad_cols)
                       if length(bad_rows) > 0 && length(bad_cols) > 0
                           outlier_cands = [bad_rows;bad_cols];
                           outlier_ind = mode(outlier_cands);
                           fprintf("dropping bad difference %d ...\n",differences(outlier_ind))
                           differences(outlier_ind) = [];
                           break_inds_neuronavi_proper(outlier_ind) = [];
                       else
                           break
                       end
                   else
                       break
                   end
               end
               offset_mean = mean(differences)
               fprintf('after checking bad differences found = %g\n',differences/1000)
           if ~all(diff(break_inds_neuronavi_proper)>0) || ~all(diff(break_starts_inds_neuronavi)>0)
                disp("break_inds_neuronavi_proper or break_starts_inds_neuronavi are not in increasing order")
                return
           end
           end
           for ind_neuronavi=1:length(break_starts_inds_neuronavi)
               if sum(ismember(break_inds_neuronavi_proper,break_starts_inds_neuronavi(ind_neuronavi)))==0
                   %not a good break index, use mean here instead
                   fprintf("setting mean offset with neuronavi index %d\n",break_starts_inds_neuronavi(ind_neuronavi))
                   offset = offset_mean;
               else
                   index_of_offset = find(break_inds_neuronavi_proper==break_starts_inds_neuronavi(ind_neuronavi));
                   offset = differences(index_of_offset);
               end
               break_inds_neuronavi(end+1) = break_starts_inds_neuronavi(ind_neuronavi);
               offsets(end+1) = offset;
           end
           string(offsets/1000)
        else
            disp("not valid neuronavi system input")
            return
        end
       %% match triggers of the stimulation/neuronavigation system with the triggers
       found_indices_all = [];
       trigger_matched_urevents_all = [];
       trigtimes = [];
       min_differences_trials_all = [];
       all_matched_trigtimes = [];

       for i=1:length(trigger_info.trigger_markers)
           %determine the current offset...
           if found_odd_offset
               offset = offsets(1);
           else
               indices_neuronavi_block_starts = find(break_inds_neuronavi(break_inds_neuronavi <= i));
               if sum(indices_neuronavi_block_starts) == 0
                   %then pick the first ever offset determined by the end of the first block
                   offset = offsets(1);
               else
                   [index_val, index_to_take] = max(indices_neuronavi_block_starts);
                   offset = offsets(index_to_take);
               end
           end
           trigger_time_in_real_world = double(trigger_info.trigger_markers{i}.recording_time_in_real_world) - offset;
           trigtimes(end+1) = trigger_time_in_real_world;
           [min_diff, argmin] = min(abs(urevents_times.times - trigger_time_in_real_world));
           trigger_urevent = urevents_times.urevents(argmin);

           if min_diff < min_difference_threshold
               if sum(ismember(trigger_matched_urevents_all,trigger_urevent)) > 0
                   disp("Matched at least one event wrongly. Checking differences in times.")
                   length(trigger_matched_urevents_all)
                   index = find(trigger_matched_urevents_all==trigger_urevent);
                   %trigger_urevent
                   %min_diff
                   %trigger_info.trigger_markers{i}.stimulus_id
                   %min_differences_trials_all(index)
                   if min_diff < min_differences_trials_all(index)
                       trigger_matched_urevents_all(index) = [];
                       min_differences_trials_all(index) = [];
                       found_indices_all(index) = [];
                       all_matched_trigtimes(index) = [];
                   end
               end
               trigger_matched_urevents_all(end+1) = trigger_urevent;
               min_differences_trials_all(end+1) = min_diff;
               found_indices_all(end+1) = i;
               all_matched_trigtimes(end+1) = trigger_time_in_real_world;
       else
           fprintf("did not find %s\n",string(trigger_info.trigger_markers{i}.stimulus_id))
        end
           

       end
     trigger_markers_remaining = trigger_info;

     if ~all(diff(trigtimes)>0) || ~all(diff(trigger_matched_urevents_all)>0) || ~all(diff(found_indices_all)>0) || ~all(diff(all_matched_trigtimes)>0)
        disp("events/indices/times are not in increasing order")
        return
     end
     fprintf("min diff and max diff %f %f\n",min(min_differences_trials_all),max(min_differences_trials_all))
     %create a dummy struct for inserting to positions that do not have
     %successful targets after alignment
     falsestruct = {struct('coil_normal',NaN,'coil_dir',NaN,'coil_ori',NaN,'coil_pos',NaN,'recording_time_latency',NaN, ...
         'recording_time_in_real_world',NaN,'matrix',NaN,'stimulus_id',NaN,'sequence_started',NaN)};

     %display how successful the matching was
     trigger_markers_remaining.trigger_markers = [];
     if length(found_indices_all) == length(trigger_info.trigger_markers)
         disp("all triggers have been matched with events in eeg all trials")
     else
         disp("all triggers have NOT been matched with events in eeg all trials")
     end
     if length(found_indices_all) == EEG_merged.trials
         disp("all trials have been matched with trigger markers")
     else
         disp("all trials have NOT been matched with trigger markers")
     end

     % fill in missing values with NaNs
     n_trials_not_matched_good = 0;
     for i=1:length(urevents_times.urevents)
         urevent_now = urevents_times.urevents(i);
         if sum(ismember(urevents_remaining_info.urevents,urevent_now))==1
             if sum(ismember(trigger_matched_urevents_all,urevent_now))==0
                trigger_markers_remaining.trigger_markers = [trigger_markers_remaining.trigger_markers,falsestruct];
                fprintf("adding falsestruct for event %d\n",urevent_now)
                n_trials_not_matched_good = n_trials_not_matched_good + 1;
             else
                index = find(trigger_matched_urevents_all==urevent_now);
                index_trigger = found_indices_all(index);
                trigger_markers_remaining.trigger_markers = [trigger_markers_remaining.trigger_markers,trigger_info.trigger_markers(index_trigger)];
             end
         elseif sum(ismember(urevents_remaining_info.urevents,urevent_now)) > 1
                disp("duplicate events, failed")
                return
         end
     end
     fprintf("eeg good trials has %d unmatched events with coil %d\n",n_trials_not_matched_good)

     %check that the times are in increasing order
     recording_times = [];
     for i =1:length(trigger_markers_remaining.trigger_markers)
         if ~isnan(trigger_markers_remaining.trigger_markers{i}.stimulus_id) %only take valid stimuli
            recording_times(end+1) = trigger_markers_remaining.trigger_markers{i}.recording_time_in_real_world;
         end
     end
     if ~all(diff(recording_times)>0)
        disp("recording_times are not in increasing order")
        return
     end

     fprintf(subject_name)
     ntrials = EEG_preprocessed.trials
     foundinds_all = length(found_indices_all)
     trigmarks_inc_nanmarkers = length(trigger_markers_remaining.trigger_markers)
     remaining_events = length(urevents_remaining_info.urevents)
     all_trigmarks = length(trigger_info.trigger_markers)

     if length(urevents_remaining_info.times) ~= length(trigger_markers_remaining.trigger_markers)
         disp("Triggers do not have trigger markers!")
         length(trigger_markers_remaining.trigger_markers)
         length(urevents_remaining_info.times)
         plot(trigtimes)
         hold on
         plot(urevents_remaining_info.times)
         return
     end
     save(trigger_markers_remaining_path,'trigger_markers_remaining')
     disp("processed")
     %return
    end
end