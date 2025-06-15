#predictions
import os
import scipy
import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from statsmodels.regression.mixed_linear_model import MixedLM
from sklearn.metrics import r2_score
from reftep_util_funcs import *
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

    

def load_data_subject(subject, site, feature_path_site, where_now, spatial_names, freq_range_names, usepsd, usecoil, usephase, usepac, phasefreqs, usetime, distance_thresh, angle_distance_thresh, angle_diff_default_normal, angle_diff_default_dir, pos_diff_default, combine_normal_ori, combine_all, usecoil_sigma):
    feature_path_subject = os.path.join(feature_path_site,subject,where_now)
    if not usepsd:
        return False #must use psd
    if usepsd:
        psds_in = {name:{} for name in spatial_names}
        for freq_range_name in freq_range_names: #go through freq_ranges
            psd_filepath_now = os.path.join(feature_path_subject,f'{subject}_{freq_range_name}_bandpowers.npy')
            psds = read_psds(psd_filepath_now)
            for name_ind, name in enumerate(spatial_names):
                psds_in[name][freq_range_name] = np.log(psds[:,name_ind]) #psds in labels/channels in trials, log-transformed
            n_trials_psd = psds.shape[0] #number of trials

    if usetime:
        if usetime[3] == "no_preprocessed_index":
            stimulation_times_path = os.path.join(feature_path_site,subject,f'{subject}_neurone_stimulation_times.mat')
            times = read_times(stimulation_times_path,usetime[1],n_trials_psd)
            if usetime[2][site] is not False and usetime[1]=='sample':
                times = times * usetime[2][site]
            shifted_times = np.log(times - np.min(times) + 1) #shift the times to account for latencies, add 1 for log transforming
            if len(shifted_times)!=n_trials_psd:
                return False
        elif usetime[3] == "preprocessed_index":
            times = np.arange(1, n_trials_psd+1, 1)
            times = times / np.max(times) #scale by the maximum
            shifted_times = np.log(times)
        else:
            return False


    if usecoil and usepsd:
        coil_control_path = os.path.join(feature_path_site,subject,f'{subject}_stimulations_final.mat')
        #read coil control params and cube-root transform them
        diff_pos, diff_normal, diff_ori, good_indices_distances = read_coil_control(coil_control_path,n_trials_psd,pos_diff_default,angle_diff_default_normal,angle_diff_default_dir,distance_thresh,angle_distance_thresh, subject)
        if min(diff_pos[good_indices_distances]) > distance_thresh or min(diff_normal[good_indices_distances]) > angle_distance_thresh or min(diff_ori[good_indices_distances]) > angle_distance_thresh:
            raise ValueError("bad distances still present in coil pos")
        if len(diff_pos) != len(good_indices_distances):
            print(f"Not using {len(diff_pos) - len(good_indices_distances)} trials due to position deviation for {subject}")
            print(f"Good coil trials {len(good_indices_distances)} out of {len(diff_pos)} for {subject}")
        #print(f'{subject} has {len(good_indices_distances)} trials in use')
        if combine_normal_ori:  #cube-root transform all variables with or without combining diff_ori and diff_normal
            diff_pos, diff_ori = np.log(diff_pos + usecoil_sigma), np.log(diff_ori + diff_normal + usecoil_sigma)
            coil_control_params = {'diff_pos':diff_pos,'diff_ori':diff_ori}
        elif combine_all: #then combine all measures after scaling each array by their maximum value
            data_coil = np.array([diff_pos, diff_normal, diff_ori])
            combined_values = np.average(data_coil, weights=np.array([0.7,0.15,0.15]), axis=0)
            coil_control_params = {'diff_coil':np.log(combined_values + usecoil_sigma)}
        else:
            diff_pos, diff_normal, diff_ori = np.log(diff_pos+ usecoil_sigma), np.log(diff_normal+ usecoil_sigma), np.log(diff_ori+ usecoil_sigma) #cube-root transform all variables
            coil_control_params = {'diff_pos':diff_pos, 'diff_normal':diff_normal,'diff_ori':diff_ori}
        if usecoil and usepsd:
            if n_trials_psd != len(diff_pos) != len(diff_ori):
                return False
    elif not usecoil and usepsd:
        good_indices_distances = np.arange(n_trials_psd)
    else:
        raise ValueError("psd must be used to get n_trials here")
    if usephase:
        for freq_range_name in phasefreqs: #go through freq_ranges
            phase_filepath_now = os.path.join(feature_path_subject,f'{subject}_{freq_range_name}_phases.mat')
            phases = read_phases(phase_filepath_now)
    if usepac:
        pacs_in = {name:{} for name in spatial_names}
        for freq_ind1, freq_range_name1 in enumerate(freq_range_names):
            for freq_ind2, freq_range_name2 in enumerate(freq_range_names):
                if freq_ind1 < freq_ind2: #only take lower-upper
                    pac_filepath_now = os.path.join(feature_path_subject,f'{subject}_{freq_range_name1}-{freq_range_name2}_pac.npy')
                    pacs = read_pacs(pac_filepath_now)
                    for name_ind, name in enumerate(spatial_names):
                        pacs_in[name][f'{freq_range_name1}_{freq_range_name2}'] = np.log(np.max(pacs[:,name_ind,:,:],axis=(1,2))) #pacs in labels/channels in trials, log-transformed
                        #print(pacs_in[name][f'{freq_range_name1}_{freq_range_name2}'].shape)


    #store the data to a coherent dictionary structure, also restrict the data to only use the good indices on distances
    subject_data = {name:{} for name in spatial_names}

    #first check for nan values in phases and remove them
    good_phase_inds = []
    bad_phase_inds = []
    for name_ind, name in enumerate(spatial_names):
        if usepsd: #save psd data
            for freq_ind, freq_range_name in enumerate(freq_range_names):
                if usephase and freq_range_name in phasefreqs:
                    for ind, (sin_phase, cos_phase) in enumerate(zip(phases[name_ind]['sines'],phases[name_ind]['cosines'])):
                        if str(sin_phase) == "nan" or str(cos_phase) == "nan":
                            bad_phase_inds.append(ind)
                        else:
                            good_phase_inds.append(ind)
    if usephase:
        good_phase_inds_unique = np.unique(good_phase_inds)
        unique_bad_phase_inds = np.unique(bad_phase_inds)
        if len(unique_bad_phase_inds) > 0:
            print(f'Found {len(unique_bad_phase_inds)} nan phases (trial indices where at least one was bad) for {subject}')
        good_phase_inds_all = np.array([ind for ind in good_phase_inds_unique if ind not in unique_bad_phase_inds])
    else:
        good_phase_inds_all = np.arange(n_trials_psd)
    good_indices_all = np.intersect1d(good_indices_distances, good_phase_inds_all) #only take valid trials
    for name_ind, name in enumerate(spatial_names):
        if usepsd: #save psd data
            for freq_ind, freq_range_name in enumerate(freq_range_names):
                subject_data[name][freq_range_name + '_PSD'] = psds_in[name][freq_range_name][good_indices_all]
                if usephase and freq_range_name in phasefreqs:
                    phases_sines = []
                    phases_cosines = []
                    for sin_phase, cos_phase in zip(phases[name_ind]['sines'],phases[name_ind]['cosines']):
                        phases_sines.append(sin_phase)
                        phases_cosines.append(cos_phase)
                    #transform to np.arrays to allow indexing as follows
                    phases_sines = np.array(phases_sines)
                    phases_cosines = np.array(phases_cosines)
                    if len(phases_sines) != len(phases_cosines) != n_trials_psd:
                        return False
                    subject_data[name][freq_range_name + '_phase_sin'] = phases_sines[good_indices_all]
                    subject_data[name][freq_range_name + '_phase_cos'] = phases_cosines[good_indices_all]
                if usepac: #save pac data to the structure
                    for freq_ind2, freq_range_name2 in enumerate(freq_range_names):
                        if freq_ind < freq_ind2: 
                            subject_data[name][f'pac_{freq_range_name}_{freq_range_name2}'] = pacs_in[name][f'{freq_range_name}_{freq_range_name2}'][good_indices_all]
        if usecoil: #these are same for all channels/labels
            for coil_param in list(coil_control_params.keys()):
                subject_data[name][coil_param] = coil_control_params[coil_param][good_indices_all]
        if usetime[0]:
            subject_data[name]['Latency'] = shifted_times[good_indices_all] #add time
        subject_data[name]['Site'] = site #site of the measurement
    #print(f"{subject} has in total {len(good_indices_all)} trials in use out of {n_trials_psd} trials")
    return subject_data



def load_responses_subject(sourcepath, subject, response):
    if response == 'mep': #get mep data and return np.inf (for fake comparisons with gof)
        amplitudes = np.log(np.load(os.path.join(sourcepath,subject,f'{subject}_{response}_amplitudes.npy')))
        return amplitudes, np.inf, np.inf
    else:
        #read all dipole amplitudes from all trials
        dipoles_trials_path = os.path.join(sourcepath, subject,f'{subject}_dipoles',f'{subject}_dipole_{response}_trials')
        n_files = len(os.listdir(dipoles_trials_path))
        amplitudes = np.array([mne.read_dipole(os.path.join(dipoles_trials_path,f'{subject}_dipole_{response}_trial_{ind}'),verbose=False).amplitude[0] for ind in range(n_files)])
        average_dipole = mne.read_dipole(os.path.join(sourcepath, subject,f'{subject}_dipoles',f'{subject}_dipole_{response}'),verbose=False) #get the gof of the dipole fitted to the average response
        average_gof = average_dipole.gof[0]
        average_amp = average_dipole.amplitude[0]
        #amplitudes = amplitudes + np.abs(np.min(amplitudes)) + 1e-9
        return np.log(np.abs(amplitudes)), average_gof, average_amp

def load_dipole_stats_subject(sourcepath, subject, response):
    #read all dipole amplitudes from all trials
    dipoles_trials_path = os.path.join(sourcepath, subject,f'{subject}_dipoles',f'{subject}_dipole_{response}_trials')
    n_files = len(os.listdir(dipoles_trials_path)) #number of files in the path
    mean_dipole = mne.read_dipole(os.path.join(sourcepath, subject,f'{subject}_dipoles',f'{subject}_dipole_{response}'),verbose=False) #the dipole fitted to the average response
    all_dipoles = [mne.read_dipole(os.path.join(dipoles_trials_path,f'{subject}_dipole_{response}_trial_{ind}'),verbose=False) for ind in range(n_files)] #dipoles fitted to single trials
    #single trial fitting information
    gofs = [dipole.gof[0] for dipole in all_dipoles]
    amplitudes = [dipole.amplitude[0] for dipole in all_dipoles]
    oris = [dipole.ori[0] for dipole in all_dipoles]
    times = [dipole.times[0] for dipole in all_dipoles]

    return gofs, mean_dipole.gof[0], amplitudes, mean_dipole.amplitude[0], oris, mean_dipole.ori[0], times, mean_dipole.times[0], mean_dipole.pos[0], mean_dipole.name, all_dipoles, mean_dipole
    


def read_coil_control(coil_control_path,n_trials,pos_diff_default,angle_diff_default_normal,angle_diff_default_dir, distance_thresh, angle_distance_thresh, subject):
    try:
        if "Aalto" in coil_control_path:
            target_name = "target_marker"
        elif 'Tuebingen' in coil_control_path:
            target_name = "instrument_marker"
        coil_control = scipy.io.loadmat(coil_control_path,simplify_cells=True)['trigger_markers_remaining']['trigger_markers']
        target = scipy.io.loadmat(coil_control_path,simplify_cells=True)['trigger_markers_remaining'][target_name]
        good_inds_notnan = np.array([i for i in range(len(coil_control)) if str(coil_control[i]['stimulus_id'])!='nan']) #some stimuli are not detected

        if 'sub-020' in coil_control_path: #here the coil control was constantly off, maybe bad instrument marker,...so using median target
            target_pos = np.median([np.array(coil_control[i]['coil_pos']) for i in good_inds_notnan],axis=0)
            target_normal = np.median([np.array(coil_control[i]['coil_normal']) for i in good_inds_notnan],axis=0)
            target_ori = np.median([np.array(coil_control[i]['coil_dir']) for i in good_inds_notnan],axis=0)
        else:
            target_pos = target['coil_pos']
            target_normal = target['coil_normal']
            target_ori = target['coil_dir']
        #get coil pos, ori and normal with dummy values for nan values
        coil_positions_all = [np.array(coil_control[i]['coil_pos']) if i in good_inds_notnan else np.array(target_pos) for i in range(len(coil_control))]
        coil_oris_all = [np.array(coil_control[i]['coil_dir']) if i in good_inds_notnan else np.array(target_ori) for i in range(len(coil_control))]
        coil_normals_all = [np.array(coil_control[i]['coil_normal']) if i in good_inds_notnan else np.array(target_normal) for i in range(len(coil_control))]

        distances_poses = np.linalg.norm(coil_positions_all - np.array(target_pos),axis=1) #distances of position to the target marker of each stim
        good_indices_distances_poses = list(np.where(distances_poses <= distance_thresh)[0])

        #angle differences to the target, get those that are within the threshold
        distances_oris = get_angle_differences(coil_oris_all, target_ori)
        good_indices_distances_oris = list(np.where(distances_oris <= angle_distance_thresh)[0])

        distances_normals = get_angle_differences(coil_normals_all, target_normal)
        good_indices_distances_normals = list(np.where(distances_normals <= angle_distance_thresh)[0])

        good_indices_distances = np.intersect1d(np.intersect1d(good_indices_distances_oris, good_indices_distances_normals),good_indices_distances_poses)

        good_inds = np.intersect1d(good_indices_distances,good_inds_notnan).astype(int) #all good indices
        median_diff_pos = np.median(distances_poses[good_inds_notnan])
        median_diff_normal = np.median(distances_normals[good_inds_notnan])
        median_diff_ori = np.median(distances_oris[good_inds_notnan])
        #replace bads with median
        diff_pos = np.array([dist if ind in good_inds else median_diff_pos for ind, dist in enumerate(distances_poses)])
        diff_normal = np.array([dist if ind in good_inds else median_diff_normal for ind, dist in enumerate(distances_normals)])
        diff_ori = np.array([dist if ind in good_inds else median_diff_ori for ind, dist in enumerate(distances_oris)])
        if subject[4] in ['1']: #all positions are recorded at aalto so if there is a nan then the pos was not available in neuronavigation
            if len(diff_pos) != len(good_inds_notnan):
                print(f"{len(diff_pos)-len(good_inds_notnan)} nan stimuli found for {subject} that is not from Tuebingen, these are rejected")
            good_indices_distances = good_inds
    except FileNotFoundError: #then replace data with the default values
        good_inds_notnan = None
        good_inds = None
        print(f"Using default difference values of position diff: {pos_diff_default} mm and normal and ori: {angle_diff_default_normal} and {angle_diff_default_dir} degrees as {coil_control_path} was not found for {subject}.")
        diff_pos = np.ones(shape=(1,n_trials))[0]*pos_diff_default
        diff_normal = np.ones(shape=(1,n_trials))[0]*angle_diff_default_normal
        diff_ori = np.ones(shape=(1,n_trials))[0]*angle_diff_default_dir
        good_indices_distances = np.arange(n_trials)
    return diff_pos, diff_normal, diff_ori, good_indices_distances #return the coil control parameters


def get_angle_differences(directions, target_direction):
    #get the angle differences of each stimulation relative to the target
    angle_list = []
    for compare in directions:
        dotprod = np.dot(compare,target_direction)
        magnitude_ori = np.linalg.norm(compare)
        magnitude_target_ori = np.linalg.norm(target_direction)
        angle = np.degrees(np.arccos(np.clip(dotprod/(magnitude_ori*magnitude_target_ori),-1.0,1.0)))
        angle_list.append(angle)
    return np.array(angle_list)

def read_psds(psd_path):
    psds_here = np.load(psd_path, allow_pickle=True)
    return psds_here

def read_phases(phase_path):
    phases = scipy.io.loadmat(phase_path, simplify_cells=True)['subject_phases']['names']
    return phases

def read_pacs(pac_path):
    pacs = np.load(pac_path)
    return pacs

def read_times(time_path, type_of_time,n_trials):
    #read the times as latencies or samples
    if type_of_time == "latency":
        times = scipy.io.loadmat(time_path,simplify_cells=True)['urevents_remaining_info']['times']
    elif type_of_time == "sample":
        #times = np.arange(n_trials)
        times = scipy.io.loadmat(time_path,simplify_cells=True)['urevents_remaining_info']['epochs_inds']
        if np.min(times) < 1 or np.max(times) > 1200:
            print("bad indices")
    else:
        raise ValueError(f'bad usetime variabe {type_of_time}, should be latency or sample')

    return np.array(times)

def scale_subject_var(group, cols_to_scale):
    scaler = StandardScaler()
    group = group.copy()
    group[cols_to_scale] = scaler.fit_transform(group[cols_to_scale])
    return group

def create_scaled_df(df, scale_within_subs, usephase):
    #scale predictive variables to zero mean and unit variance
    scaler = StandardScaler()
    #pick only some explanatory variables
    explanatory_vars_to_scale = [col for col in df.columns if col not in ['AMP', 'Sample', 'Subject', 'Site', 'phase_sin_alpha', 'phase_cos_alpha'] and col not in scale_within_subs]
    df_scaled = df.copy()
    #fit the transform globally to some variables
    df_scaled[explanatory_vars_to_scale] = scaler.fit_transform(df[explanatory_vars_to_scale])
    #scale subject-wise to some variables
    df_scaled = df_scaled.groupby('Subject').apply(lambda g: scale_subject_var(g, scale_within_subs)).reset_index(drop=True)
    numeric_vars_to_use = explanatory_vars_to_scale + scale_within_subs
    if usephase:
        numeric_vars_to_use += ['phase_sin_alpha', 'phase_cos_alpha']
    return df_scaled, numeric_vars_to_use



def create_df_with_features(data, responses, name, freq_range_names, phases, pac, coil_control, phasefreqs, usetime, combine_normal_ori, combine_all):
    n_subjects = len(data) #number of subjects
    list_for_data = []
    for subject_idx in range(n_subjects):
        subject_name = f'S{subject_idx+1}'
        all_keys =list(data[subject_idx][name].keys())
        first_key = all_keys[0]
        if len(np.unique([len(data[subject_idx][name][key]) for key in all_keys if key != 'Site'])) != 1:
            print("failed")
        n_samples = len(data[subject_idx][name][first_key]) #number of samples
        site = data[subject_idx][name]['Site']
        for sample_idx in range(n_samples): #go through all samples
            row = {'Subject': subject_name, 'Sample': sample_idx, 'AMP':responses[subject_idx][sample_idx], 'Site':site} #initialize row with current subject and the sample
            if freq_range_names: #include frequency ranges (include power)
                for freq_ind, freq_range_name in enumerate(freq_range_names): #add info from PSDs
                    row[f'PSD_{freq_range_name}'] = data[subject_idx][name][f'{freq_range_name}_PSD'][sample_idx]
                    """if connectivity: #add connectivity data
                        for pos_ind, (name2, connectivity) in enumerate(data[subject_idx][name]['cons'][freq_range_name].items()):
                            row[f'CON_{freq_range_name}_{name2}'] =  connectivity[sample_idx]
                        return False"""
                    if phases and freq_range_name in phasefreqs: #include estimated signal phases at the time of the tms pulse
                        row[f'phase_sin_{freq_range_name}'] = data[subject_idx][name][f'{freq_range_name}_phase_sin'][sample_idx]
                        row[f'phase_cos_{freq_range_name}'] = data[subject_idx][name][f'{freq_range_name}_phase_cos'][sample_idx]
                    if pac: #include phase-amplitude coupling
                        for freq_ind2, freq_range_name2 in enumerate(freq_range_names):
                            if freq_ind < freq_ind2: #only take lower-upper
                                row[f'pac_{freq_range_name}_{freq_range_name2}'] = data[subject_idx][name][f'pac_{freq_range_name}_{freq_range_name2}'][sample_idx]
                                row[f'pac_{freq_range_name}_{freq_range_name2}'] = data[subject_idx][name][f'pac_{freq_range_name}_{freq_range_name2}'][sample_idx]
            if coil_control: #include coil control paramaters
                if combine_normal_ori:
                    row['diff_pos'] = data[subject_idx][name]['diff_pos'][sample_idx]
                    row['diff_ori'] = data[subject_idx][name]['diff_ori'][sample_idx]
                elif combine_all:
                    row['diff_coil'] = data[subject_idx][name]['diff_coil'][sample_idx]
                else:
                    row['diff_pos'] = data[subject_idx][name]['diff_pos'][sample_idx]
                    row['diff_normal'] = data[subject_idx][name]['diff_normal'][sample_idx]
                    row['diff_ori'] = data[subject_idx][name]['diff_ori'][sample_idx]
            if usetime[0]:
                row['Latency'] = data[subject_idx][name]['Latency'][sample_idx]
            list_for_data.append(row)
    
    df = pd.DataFrame(list_for_data)
    #for freq_range_name in freq_range_names: #baseline-correct for each subject individually
        #df[f'PSD_{freq_range_name}_mean'] = df.groupby('Subject')[f'PSD_{freq_range_name}'].transform('mean')
        #df[f'PSD_{freq_range_name}'] = df[f'PSD_{freq_range_name}'] - df[f'PSD_{freq_range_name}_mean']
        #del df[f'PSD_{freq_range_name}_mean'] #don't use mean
    return df



def perform_lmm(df, explanatory_variables, grouptype, interaction_variables_with_time, ref_site, re_formula=False):
    if ref_site is not None:
        formula = 'AMP ~ C(Site) + ' + ' + '.join(explanatory_variables) #AMP ~ Feature_1 + Feature_2 +... + Feature_n for n features
    else:
        formula = 'AMP ~ ' + ' + '.join(explanatory_variables) #AMP ~ Feature_1 + Feature_2 +... + Feature_n for n features
    #add interaction effects
    interaction_vars = []
    for var in interaction_variables_with_time:
        interaction_vars += [var + '_x_Latency']
        df[var + '_x_Latency'] = df[var] * df['Latency']
    if len(interaction_variables_with_time) > 0:
        formula = formula + '+ ' + ' + '.join(interaction_vars)
    df_new = df.copy()
    if ref_site is not None:
        new_ref_site = "AAAA" + ref_site
        df_new['Site'] = df_new['Site'].replace(ref_site,new_ref_site)
    if re_formula:
        model = MixedLM.from_formula(formula, df_new, groups=grouptype, re_formula=re_formula) #create the model with random effect
    else:
        model = MixedLM.from_formula(formula, df_new, groups=grouptype) #create the model with random effect
    result = model.fit() #fit the model
    return result



def calculate_r2(df, target_name, result):
    #calculate explained variance
    return r2_score(df[target_name], result.fittedvalues)

def calculate_r2_cond_marginal(result, df):
     #calculate marginal and conditional r-squared
    #var_fitted = np.var(result.fittedvalues) #variance of fitted values
    exog = result.model.exog
    fe_params = result.fe_params
    sigma_f2 = np.var(np.dot(exog, fe_params.values),ddof=1) #variance of fixed effects

    #fixed_effects_pred = result.predict(exog=df)
    #sigma_f2 = np.var(fixed_effects_pred) #variance of fixed effects

    sigma_r2 = result.cov_re.iloc[0,0] #variance of random effects
    sigma_e2 = result.scale #residual error variance
    #print(sigma_f2, sigma_r2, sigma_e2)
    r2_conditional = (sigma_f2 + sigma_r2) / (sigma_f2 + sigma_r2 + sigma_e2)
    r2_marginal = sigma_f2 / (sigma_f2 + sigma_r2 + sigma_e2)
    return r2_conditional, r2_marginal


def estimate_phase_coef(sin_param, cos_param, model):
    b_sin = model.params[sin_param]
    b_cos = model.params[cos_param]
    #use sin and cos coefs to estimate coef of phase
    b_phase = np.sqrt(b_sin**2 + b_cos**2)*np.sign(b_sin)
    if b_phase == 0:
        print("b_phase is 0")
        return b_phase, np.nan
    vcovs = model.cov_params()
    #extract variances of sin and cos
    var_bsin = vcovs.loc[sin_param, sin_param]
    var_bcos = vcovs.loc[cos_param, cos_param]
    #extract covariance
    var_bsin_bcos = vcovs.loc[sin_param, cos_param]
    #approximate the variance
    var_bphase = var_bsin*((b_sin**2)/(b_phase**2)) + var_bcos*((b_cos**2)/(b_phase**2)) + 2*var_bsin_bcos*((b_sin*b_cos)/(b_phase**2))
    se_phase = np.sqrt(var_bphase) #standard error
    z_score_res = b_phase/se_phase #z-score, Wald test
    p_phase = 2*(1-norm.cdf(np.abs(z_score_res)))
    return b_phase, p_phase

def perform_p_val_correction(result, p_phase):
    p_values_nointercept = result.pvalues.drop("Intercept")
    if p_phase is not False:
        p_values = pd.Series(np.append(p_values_nointercept.values,p_phase), index=p_values_nointercept.index.append(pd.Index(['phase_alpha']))) #don't use the intercept but add phase
    else:
        p_values = p_values_nointercept
    #adjust for multiple comparisons
    corrected_pvals_bonferroni = multipletests(p_values.values, method='bonferroni')[1]
    corrected_pvals_fdr = multipletests(p_values.values, method='fdr_bh')[1]
    corrected_p_values = pd.DataFrame({'predictor': p_values.index,'raw_p_val':p_values, 'bonferroni':corrected_pvals_bonferroni,'fdr':corrected_pvals_fdr})
    return corrected_p_values










