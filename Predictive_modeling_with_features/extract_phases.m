%estimate phases
addpath("phastimate-master\"); %phastimate
addpath("eeglab2024.2\");
eeglab nogui;

fieldnames_freqs = {'alpha'};
spatial = {"['n15', 'p30', 'n45', 'p60', 'handknob']",'aparc'};

ar_order = 25;
edge = 65;
hilbert_window = 128;
offset_correction = 46;
filter_order = 192;
sampling_rate = 1000;
%n_samples = [1000];

%what to plot?
plot_distribution = false;
plot_name = "['n15', 'p30', 'n45', 'p60', 'handknob']";
plot_index = 5;

for spatial_ind = 1:length(spatial)
    where = char(spatial{spatial_ind}); %sensors or sources
    for site={'Aalto','Tuebingen'}
        site_char = char(site);
        datapath = char(strcat('D:\REFTEP_ALL\Source_analysis\Source_analysis_',site_char));
        datapath_features = char(strcat('D:\REFTEP_ALL\Features\Features_',site_char,'\'));
        files_and_folders = dir(datapath);
        is_subfolder = [files_and_folders.isdir];
        folders = files_and_folders(is_subfolder);
        names = {folders.name};
        subject_names = names(contains(names,"sub"));
        for index = 1:length(subject_names)%go through all subjects
            subject = subject_names(index); %name of the subject
            directory_name = fullfile(datapath_features,subject); %features of the subject in this path
            %load the individual freq_ranges based on the alpha peak
            freq_range_filename = char(strcat(subject,'-freq_ranges_dict_matlab.mat'));
            freq_ranges = load(char(fullfile(directory_name,freq_range_filename)));
            for freq_ind = 1:length(fieldnames_freqs)
                freq_range_now = freq_ranges.(fieldnames_freqs{freq_ind}); %used freq range now
                %load the parcellated source estimates
                stc_name = strcat(subject,'_',where,'_cropped_mne.mat'); %not really cropped even though the name says so
                stc_filename = char(fullfile(datapath,subject,stc_name));
                full_data = load(stc_filename).source_estimate;
                data = permute(full_data,[2,3,1]); %rearrange to n_labels x n_times x n_trials
                number_of_names = size(data,1); %number of labels
                if size(data,2) ~= 1000
                    disp("not 1000 samples, beware...")
                end
                %data_cropped = data(:,end-n_samples(freq_ind)+1:end,:); %pick n samples
                D_bandpass = designfilt('bandpassfir', 'FilterOrder', filter_order, 'CutoffFrequency1', freq_range_now(1) , 'CutoffFrequency2', freq_range_now(2), 'SampleRate', sampling_rate, 'DesignMethod', 'window');
                subject_phases = struct();
                for index_of_loc=1:number_of_names
                    signal = squeeze(data(index_of_loc,:,:)); %squeeze the data of one label
                    %estimate the phases of the cycle at the stimulus
                    [phases, amplitudes] = phastimate(signal, D_bandpass, edge, ar_order, hilbert_window, offset_correction);
                    if index_of_loc == plot_index && strcmp(where,plot_name) && plot_distribution
                        phases_deg = rad2deg(phases); %ransform to degrees
                        figure; %plot figure
                        polarhistogram(phases_deg,30)
                        title(strcat(fieldnames_freqs{freq_ind},'-',where,'-',char(subject)))
                        break
                    end
                    %get the sin and cos of the phase
                    sines = sin(phases);
                    cosines = cos(phases);
                    %create a sub structure with the extracted information
                    subject_phases.names{index_of_loc} = struct('sines', sines, 'cosines', cosines, 'phases', phases);
                end
                %save the struct for the phases for each channel
                subject_phases_path = char(fullfile(directory_name,char(strcat('source/',subject,'_',where,'/',subject, '_', fieldnames_freqs{freq_ind},'_phases.mat'))));
                save(subject_phases_path, "subject_phases"); %save the phases and the related info
            end
        end
    end
end