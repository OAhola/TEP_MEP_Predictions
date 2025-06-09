function epoch_inds = get_epoch_inds_good(EEG,exclude_range)
    epoch_inds = [];
    n_epochs = size(EEG.data,3);
    if n_epochs ~= size(EEG.epoch,2)
        disp("failed size")
        return
    end
    for i=1:n_epochs
        epoch_event = EEG.epoch(i).event;
        if epoch_event ~= i
            disp("failed index")
            epoch_inds = [];
            return
        end
        if ~ismember(i, exclude_range)
            epoch_inds(end+1) = i;
        else
            i
            exclude_range
        end
    end
end