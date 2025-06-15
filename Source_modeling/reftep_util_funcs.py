import mne
import numpy as np
import sklearn.metrics
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


#Functions for dipole fitting
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def r2_dipole(forward, evoked, tmin=None, tmax=None):
    """
    forward: the forward solution
    evoked: evoked structure
    tmin, tmax = minimum and maximum time, can be tmin=tmax
    """
    L = forward['sol']['data']
    n_dips_x_oris = int(L.shape[1]) #number of dipoles x orientations
    L = L - np.mean(L,axis=0) #average reference the leadfield
    r2score = -np.inf #initialize that the r2 is -np.inf
    best_topo = None
    best_free_ori_stc_now = None
    best_time = None
    best_pos_ind = None
    best_match = None
    evoked_full_data = evoked.data
    n_channels = evoked.data.shape[0] #number of channels
    #if tmin and tmax have not been specified, then use the whole time range of the evoked response (or single trial)
    if tmin is None and tmax is None:
        evoked_cropped = evoked
    elif (tmin is not None and tmax is None) or (tmin is None and tmax is not None):
        print("tmin and tmax must both be None or not None")
        return False
    else:
        evoked_cropped = evoked.copy().crop(tmin, tmax)
    for time_ind, time in enumerate(evoked_cropped.times): #go through the times
        #y_measured = evoked_cropped.data[:,time_ind] - np.mean(evoked_cropped.data[:,time_ind], axis=0) #the measured topography at the current time in average reference
        #data is alrady in ave ref
        y_measured = evoked_cropped.data[:,time_ind]
        for lf_ind, dipole_start_index in enumerate(np.arange(0,int(n_dips_x_oris),3)): #go through lead field triplets, i.e. positions with orientations
            pickcol = L[:,dipole_start_index:dipole_start_index+3]#picked triplet here
            free_orientation_stc_now = np.matmul(np.linalg.pinv(pickcol),y_measured) #q_hat = pinv(L)y
            y_predicted = np.matmul(pickcol,free_orientation_stc_now) #y_hat = Lq_hat
            r2_new = sklearn.metrics.r2_score(y_measured, y_predicted) #r squared between y and y_hat
            if r2_new > r2score: #if the best score was exceeded, then save the dipole info as the best one
                r2score = r2_new #best new score
                #update best position
                best_pos_ind = lf_ind
                best_match = dipole_start_index
                best_topo = y_predicted #save the best (so far) topography
                best_free_ori_stc_now = free_orientation_stc_now #dipole orientation
                best_time = time #the time of the fit


    evoked_fit = mne.EvokedArray(data=best_topo.reshape(n_channels,1), info=evoked.info, tmin=best_time, nave=1, kind='average', baseline=None) #the dipolar acitity
    #project out the dipolar activity from the data and return the resulting residual
    lf_cols = L[:,best_match:best_match+3] #best matching source
    tc_of_source = np.matmul(np.linalg.pinv(lf_cols),evoked_full_data) #full time course of that source
    source_projected = np.matmul(lf_cols,tc_of_source) #time course of the source projected to sensor space
    full_residual_data = evoked_full_data - source_projected #calculate the residual (dipolar activity removed)
    #full_residual_data = full_residual_data - np.mean(full_residual_data,axis=0) #average reference, not applied here
    evoked_residual = mne.EvokedArray(data=full_residual_data, info=evoked.info, tmin=evoked.times[0], nave=1, kind='average', baseline=None) #residual data

    #extract dipole information and create a dipole object
    dipole_amplitude = np.linalg.norm(best_free_ori_stc_now) #dipole amplitude
    dipole_pos = forward['source_rr'][best_pos_ind] #dipole position in mri
    dipole_ori = best_free_ori_stc_now / dipole_amplitude #dipole orientation normalized (q/amplitude)
    dipole = mne.Dipole([best_time], [dipole_pos], [dipole_amplitude], [dipole_ori], [r2score*100], name=f"dipole_{best_pos_ind}") #create dipole object (mne.Dipole)
    return dipole, best_free_ori_stc_now, evoked_fit, best_match, best_pos_ind, evoked_residual

def dipole_to_pos(forward, evoked, tmin, tmax, pos_ind, n_times, maximize, ori_fixed=None):
    """
    This function fits a dipole to a pre-defined location and gives it a free orientation, a fixed orientation is given if ori_fixed is not None
    """
    start_of_triplet_ind = pos_ind*3
    n_channels = evoked.data.shape[0] #number of channels
    L = forward['sol']['data']
    L = L - np.mean(L,axis=0) #average reference the leadfield
    #initialize values
    r2score = -np.inf
    amplitude = 0
    best_topo = None
    best_ori_stc_now = None
    best_time = None
    evoked_cropped = evoked.copy().crop(tmin, tmax)
    if len(evoked_cropped.times) != n_times:
        print(len(evoked_cropped.times), n_times)
        return False
    pickcol = L[:,start_of_triplet_ind:start_of_triplet_ind+3] #picked triplet here
    for time_ind, time in enumerate(evoked_cropped.times):
        y_measured = evoked_cropped.data[:,time_ind]
        #the data is already average referenced
        #y_measured = evoked_cropped.data[:,time_ind] - np.mean(evoked_cropped.data[:,time_ind], axis=0) #the measured topography at the current time in average reference
        if ori_fixed is None:
            orientation_stc_now = np.matmul(np.linalg.pinv(pickcol),y_measured) #q_hat = pinv(L)y
            amplitude_now = np.linalg.norm(orientation_stc_now)
        else:
            pickcol_projected = np.matmul(pickcol,ori_fixed) #project the leadfield to the fixed orientation
            #amplitude_now = np.matmul(np.linalg.pinv(pickcol_projected),y_measured)[0] #amplitude of the fixed source, to use this you must have pickcol_projected = np.matmul(pickcol,ori_fixed)[:, np.newaxis]
            amplitude_now = np.dot(pickcol_projected.T,y_measured) /  np.dot(pickcol_projected.T, pickcol_projected) #amplitude of the fixed source
            orientation_stc_now = ori_fixed*amplitude_now #the source strengths in all dimensions
        y_predicted = np.matmul(pickcol,orientation_stc_now) #y_hat = Lq_hat
        r2_new = sklearn.metrics.r2_score(y_measured, y_predicted) #r squared between y and y_hat
        threshold_exceeded = False
        if maximize == "r2":
            if r2_new > r2score:
                threshold_exceeded = True
        else:
            raise ValueError("maximize should be r2, others are not supported currently")
        if threshold_exceeded: #if the best score (r2) was exceeded, then save the dipole info as the best one
            r2score = r2_new
            best_topo = y_predicted
            best_ori_stc_now = orientation_stc_now
            best_time = time
            amplitude = amplitude_now

    #extract dipole information and create a dipole object
    evoked_fit = mne.EvokedArray(data=best_topo.reshape(n_channels,1), info=evoked.info, tmin=best_time, nave=1, kind='single_epoch', baseline=None)
    dipole_pos = forward['source_rr'][pos_ind] #dipole position in mri
    if ori_fixed is None:
        dipole_ori = best_ori_stc_now / amplitude #dipole orientation normalized (q/amplitude)
    else:
        dipole_ori = ori_fixed
    dipole = mne.Dipole([best_time], [dipole_pos], [amplitude], [dipole_ori], [r2score*100], name=f"dipole_{pos_ind}") #create dipole object (mne.Dipole)
    return dipole, best_ori_stc_now, evoked_fit


def get_good_dipole_indices(dipole_fixeds, orig_dipole, angle_diff_thresh, gof_perc_thresh):
    good_indices = []
    for dipole_index, dipole_fixed in enumerate(dipole_fixeds):
        dotprod = np.dot(dipole_fixed.ori[0], orig_dipole.ori[0])
        #these are very near one because of normalization but due to rounding errors etc so one could divide once again (shouldn't hurt)
        #magnitude_ori_fixed = np.linalg.norm(dipole_fixed.ori[0])
        #magnitude_ori_orig = np.linalg.norm(orig_dipole.ori[0]) #they are already unit length because they have been normalized before
        #magnitude_ori_fixed = 1
        #magnitude_ori_orig = 1
        angle_diff = np.degrees(np.arccos(np.clip(dotprod,-1.0,1.0)))
        if angle_diff <= angle_diff_thresh and dipole_fixed.gof[0] >= orig_dipole.gof[0]*gof_perc_thresh:
            good_indices.append(dipole_index)
    return good_indices
    
def get_fitting_times(possible_times, best_dipole_index, good_indices):
    final_times = [] #get the times around the best dipole where the criterions match
    ind1 = 0 #start from 0 to include the best dipole time
    final_good_indices = []
    while True: #get the times at best dipole and before
        index_to_look_for_1 = best_dipole_index - ind1
        if index_to_look_for_1 in good_indices:
            final_times.append(possible_times[index_to_look_for_1])
            final_good_indices.append(index_to_look_for_1)
        else:
            break
        ind1 += 1
    ind2 = 1
    while True: #times after the best dipole
        index_to_look_for_2 = best_dipole_index + ind2
        if index_to_look_for_2 in good_indices:
            final_times.append(possible_times[index_to_look_for_2])
            final_good_indices.append(index_to_look_for_2)
        else:
            break
        ind2 += 1
    final_good_indices = np.array(sorted(final_good_indices))
    final_times = np.array(sorted(final_times)) #the final times to consider
    return final_times, final_good_indices

    
def get_mep_ptps(emg_epochs,tmin,tmax):
    data = emg_epochs.copy().crop(tmin=tmin, tmax=tmax).get_data(copy=True)
    mins = np.min(data, axis=2)
    maxs = np.max(data, axis=2)
    min_max_diffs = maxs - mins
    mean_min_max_diffs = np.mean(min_max_diffs,axis=0)
    maxind = np.argmax(mean_min_max_diffs) #which one has on average a higher amplitude
    return min_max_diffs[:,maxind], maxind #return those amplitudes (min-max differences)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




#get PSDs for defining frequency ranges based on alpha peak (if found)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def get_psd_over_freqs(epochs,fmin,fmax,bw_scale):
    sfreq = epochs.info["sfreq"] #sampling frequency
    data = epochs.get_data(copy=True) #epoch data
    n_times = data.shape[-1] #number of time points
    bw = bw_scale * (sfreq / n_times) #used bandwidth
    psds, freqs = mne.time_frequency.psd_array_multitaper(data, sfreq, fmin=fmin, fmax=fmax, bandwidth=bw, output='power')
    return psds, freqs

def get_individual_alpha_peak(epochs,bw_scale,fmin,fmax):
    psd, freqs = get_psd_over_freqs(epochs,fmin=fmin,fmax=fmax,bw_scale=bw_scale) #psd from fmin to fmax
    mean_psd = np.mean(psd, axis=0)
    mean_psd = np.mean(mean_psd, axis=0) #mean psd over epochs and channels
    peaks, _ = find_peaks(mean_psd) #find peaks in the psd
    if len(peaks) == 0:
        return False, freqs #no individual peak found
    if len(peaks) == 1:
        max_peak_ind = 0
        alpha_peak = freqs[peaks[max_peak_ind]] #peak found!
    else:
        peaks_vals = mean_psd[peaks]
        max_peak_ind = np.argmax(peaks_vals)
        alpha_peak = freqs[peaks[max_peak_ind]]
        print("more than one peak found, picked one with highest value")
    fig, axs  = plt.subplots()
    axs.plot(freqs,mean_psd)
    axs.axvline(alpha_peak, mean_psd[max_peak_ind])
    plt.show()
    return alpha_peak, freqs

def get_freq_ranges_based_on_alpha_peak(alpha_peak, freq_range_names):
    alpha = [alpha_peak-2.5, alpha_peak+2.5]
    min_alpha = np.min(alpha)
    theta = [min_alpha-3.5, min_alpha]
    max_alpha = np.max(alpha)
    beta = [max_alpha, max_alpha+20]
    max_beta = np.max(beta)
    gamma = [max_beta, max_beta+50]
    frequency_ranges = [theta, alpha, beta, gamma]
    if len(frequency_ranges) != len(freq_range_names):
        return False
    freq_range_dict = {name:freq_range for name, freq_range in zip(freq_range_names, frequency_ranges)}
    return frequency_ranges, freq_range_dict

#------------------------------------------------------------------------------------------------------------------------------------------------------------------     
