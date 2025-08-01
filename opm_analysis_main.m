%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    server = true;
    % Server:
    base_data_path = '/archive/23108_CAPSI/MEG/';
    base_save_path = '/home/share/capsi_share';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/opm_general';
else
    server = false;
    % Laptop:
    base_data_path = '/Volumes/dataarchvie/CHOP/MEG';
    base_save_path = '/Users/christophpfeiffer/data_local/CHOP';
    base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/christophpfeiffer/opm_general';
end

%% Set up fieldtrip
addpath(fullfile(base_matlab_path,'fieldtrip')) % Fieldtrip path
addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Analyse squid data?
squid = false;

%% Overwrite
overwrite = [];
overwrite.preproc = true;
overwrite.timelock = true;
overwrite.mri = false;
overwrite.coreg = false;
overwrite.dipole = false;
overwrite.mne = false;


%% Params
params = [];
params.pre = 0.3; % Trial prestim in seconds
params.post = 0.1; % Trial poststim in seconds
params.pad = 0.2; % Trial (pre and post) padding in seconds
params.delay = 0.01; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).

params.filter = [];
params.filter.hp_freq = 1; % Highpass cutoff frequency
params.filter.lp_freq = 50; % Lowpass cutoff frequency
params.filter.notch = sort([50 60 100 120]); % Notch (bandstop) filter frequencies

params.ds_freq = 1000; % Downsample frequency. If empty or not defined no downsampling will be applied

params.apply_hfc = true; % Apply Homogenous Field Correction
params.hfc_order = 1; % Order for Homogenous Field Correction: 1..3

params.apply_amm = false; % Apply Adaptive Multipole Models
params.amm_in = 12;
params.amm_out = 2;
params.amm_thr = 1;

params.n_comp = 40; % Number of ICA components 
params.manual_ica = false; % Manually select ICA components to remove?
params.save_ica = 1; % Save plots and components
params.ica_cor = 0.8; % Cutoff for correlation with EOG/ECG 
params.ica_coh = 0.95; % Cutoff for coherence with EOG/ECG 

params.corr_threshold = 0.7; % Correlation threshold for badchannel neighbors
params.z_threshold = 20; % Zmax threshold for badchannel and trial detection
params.opm_std_threshold = 5e-12; % Stddev threshold for badtrial detection
params.squid_std_threshold = 2.5e-12; % Stddev threshold for badtrial detection

params.hpi_freq = 33; % HPI coil frequency
params.hpi_gof = 0.9; % Minimum goodness-of-fit for including coil in hpi analysis

params.trigger_codes = [3 5 9]; % Trigger values to timelock
params.trigger_labels = {'short', 'middle', 'long'}; % Labels corresponding to the trigger values

params.src_density = '8'; % Sourcemodel density ('4', '8' or '32') = approximate number of sources per hemisphere
params.source_fixedori = true; % use fixed orientation sources (along vertex normals); if false: use three orthogonal sources per location
params.use_cov = 'resting_state'; % noise cov to use; default= ' ' for prestim, alt: 'resting_state', 'empty_room'

params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';

%% Subjects + dates
[subjects, sessions] = getSubjectsAndSessions(base_data_path,false);

paradigms = {'varITI'}; % Paradigms to analyze for all participants and sessions

if server
    subs_to_run = [find(cellfun(@(x) strcmp(x,'1196'), subjects)) find(cellfun(@(x) strcmp(x,'1206'), subjects)) find(cellfun(@(x) strcmp(x,'1211'), subjects))];
else
    %subs_to_run = find(cellfun(@str2num, subjects)>1082)';
    subs_to_run = find(cellfun(@(x) strcmp(x,'0953'), subjects));
end
ses_cnt = 0;

%% Loop over subjects
for i_sub = subs_to_run
    params.sub = ['sub-' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path, params.sub);  
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end
    %% Loop over sessions
    for i_ses = 1:length(sessions(i_sub,:))
        if isempty(sessions{i_sub,i_ses})
            disp(['No session defined! Skipping sub-' num2str(i_sub,'%02d') '_ses-' num2str(i_ses,'%02d')])
            continue % Skip iteration if no session defined
        end
        ses_cnt = ses_cnt + 1;

        params.ses = ['ses-' num2str(i_ses,'%02d')];
        
        %% Paths
        raw_path = fullfile(base_data_path, ['NatMEG_' subjects{i_sub}], sessions{i_sub,i_ses});
        save_path = fullfile(base_save_path, params.sub, params.ses);
        if ~exist(save_path, 'dir')
           mkdir(save_path)
        end
        if ~exist(fullfile(save_path,'figs'), 'dir')
           mkdir(fullfile(save_path,'figs'))
        end
        for i_paradigm = 1:length(paradigms)
            if server % on server
                tmp = dir(fullfile(raw_path,'opm',['*' paradigms{i_paradigm} 'OPM_raw.fif']));
                opm_files{i_paradigm} = fullfile(tmp.folder,tmp.name); % opm files 
                squid_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'MEG_tsss_mc.fif']); % corresponding aux files containing EOG/ECG
            else % on laptop
                opm_files{i_paradigm} = fullfile(raw_path,'osmeg',[paradigms{i_paradigm} 'OPM_raw.fif']); % opm files 
                squid_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']); % corresponding aux files containing EOG/ECG
            end
            aux_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'EEG.fif']); % corresponding aux files containing EOG/ECG
        end
        hpi_path = fullfile(raw_path,'osmeg');
        
        %% Loop over paradigms/tasks
        for i_paradigm = 1:length(paradigms)
            params.paradigm = paradigms{i_paradigm};
            
            %% Read and preproc
            params.modality = 'opm';
            params.layout = 'fieldlinebeta2bz_helmet.mat';
            params.chs = '*_b*';
            
            if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica.mat']),'file')
                ft_hastoolbox('mne', 1);
    
                % Read data 
                disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigms)) '...'])
                data_epo = read_osMEG(opm_files{i_paradigm}, aux_files{i_paradigm}, save_path, params); % Read data
                
                % ICA
                disp('Running ICA ...')
                if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1 % No ExG data
                    params.manual_ica = 1;
                    params.save_ica = 1;
                end
                data_ica = ica_MEG(data_epo, save_path, params);
                save(fullfile(save_path, [params.paradigm '_data_ica']), 'data_ica',"-v7.3"); disp('done');
                clear data_epo
            else
                data_ica = load(fullfile(save_path, [params.paradigm '_data_ica.mat'])).data_ica;
            end
            
            if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
                params.modality = 'opm';
                params.layout = 'fieldlinebeta2bz_helmet.mat';
                params.chs = '*bz';
                params.amp_scaler = 1e15;
                params.amp_label = 'B [fT]';
                timelocked = timelock(data_ica, save_path, params);
                save(fullfile(save_path, [params.paradigm '_timelocked']), 'timelocked', '-v7.3'); 
                clear timelocked
            end
            clear data_ica
    
            if squid
                %% Read and preproc - SQUID-MAG
                params.modality = 'squid';
                params.layout = 'neuromag306mag.lay';
                params.chs = 'meg';
        
                if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat']),'file')
                    ft_hastoolbox('mne', 1);
        
                    % Read data 
                    disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigms)) '...'])
                    data_epo = read_cvMEG(squid_files{i_paradigm}, params); % Read data
                    
                    % ICA
                    disp('Running ICA ...')
                    if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1 % No ExG data
                        params.manual_ica = 1;
                        params.save_ica = 1;
                    end
                    data_ica = ica_MEG(data_epo, save_path, params);
                    save(fullfile(save_path, [params.paradigm '_data_ica_squidmag']), 'data_ica',"-v7.3"); disp('done');
                    clear data_epo
                else
                    data_ica = load(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat'])).data_ica;
                end
                
                if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
                    params.modality = 'squidmag';
                    params.layout = 'neuromag306mag.lay';
                    params.chs = 'megmag';
                    params.amp_scaler = 1e15;
                    params.amp_label = 'B [fT]';
                    timelocked = timelock(data_ica, save_path, params);
                    save(fullfile(save_path, [params.paradigm '_timelocked_squidmag']), 'timelocked', '-v7.3'); 
                    clear timelocked
                end
                clear data_ica
            end
        end
        
        %% HPI localization
        ft_hastoolbox('mne',1);
        
        if exist(fullfile(save_path_mri, 'opm_trans.mat'),'file') && overwrite.coreg==false
            disp(['Not overwriting OPM transform for ' params.sub]);
        else
            ft_hastoolbox('mne', 1);
            params.include_chs = load(fullfile(save_path, ['include_chs' num2str(length(opm_files))])).include_chs;
            clear data_ica
            opm_trans = fit_hpi(hpi_path, aux_files{1}, save_path, params);
        
            opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked;
            for i = 1:length(params.trigger_labels)
                opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.chanpos);
                opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.coilpos);
                opm_timelockedT{i}.grad.chanori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.chanori')';
                opm_timelockedT{i}.grad.coilori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.coilori')';
            end
    
            % Plot source and head models
            clear headmodels sourcemodel
            headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
            sourcemodel = load(fullfile(save_path, 'sourcemodel.mat')).sourcemodel;
        
            h=figure; 
            ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
            hold on; 
            ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
            ft_plot_sens(opm_timelockedT.grad,'unit','cm')
            hold off;
            title('OPM-MEG')
            view([-140 10])
            saveas(h, fullfile(save_path, 'figs', 'opm_layout.jpg'))
            close all
            
            save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
            clear timelocked sourcemodel headmodels opm_trans
        end

        %% Dipole fits
        ft_hastoolbox('mne',1);  
        if exist(fullfile(save_path, [params.peaks{1}.label '_dipoles.mat']),'file') && overwrite.dip==false
            disp(['Not overwriting dipole source reconstruction for ' params.sub]);
        elseif exist(fullfile(save_path, [params.sub '_opm_timelockedT.mat']),'file')
            headmodel = load(fullfile(save_path_mri, 'headmodels.mat')).headmodels.headmodel_meg;
            mri_resliced = load(fullfile(save_path_mri, 'mri_resliced.mat')).mri_resliced;
            opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
            squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked;
            
            for i_peak = 1:length(params.peaks)
                peak_opm = load(fullfile(save_path, [params.sub '_opm_' params.peaks{i_peak}.label])).peak; 
                fit_dipoles(save_path, opm_timelockedT, headmodel, mri_resliced, params);
                peak_squid = load(fullfile(save_path, [params.sub '_squid_' params.peaks{i_peak}.label])).peak; 
                fit_dipoles(save_path, squid_timelocked, headmodel, mri_resliced, params);
                clear peak_opm peak_squid
            end
            clear squid_timelocked opm_timelockedT
        end
        
        %% MNE
        ft_hastoolbox('mne',1);
        if exist(fullfile(save_path, 'opm_mne_peaks.mat'),'file') && overwrite.mne==false
            disp(['Not overwriting MNE source reconstruction for ' params.sub]);
        elseif exist(fullfile(save_path, [params.sub '_opm_timelockedT.mat']),'file') 
            clear headmodel sourcemodel sourcemodel_inflated
            sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel'])).sourcemodel;
            sourcemodel_inflated = load(fullfile(save_path, [params.sub '_sourcemodel_inflated'])).sourcemodel_inflated;
            headmodel = load(fullfile(save_path,'headmodels.mat')).headmodels.headmodel_meg;
    
            %% OPM
            clear opm_timelockedT
            opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
            
            for i = 1:length(opm_timelockedT)
                opm_timelockedT{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_opm.mat'])).opm_RS_cov;
                if exist(fullfile(save_path, [params.sub '_ER_squid.mat']),'file')
                    opm_timelockedT{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_opm.mat'])).opm_ER_cov;
                end
            end
            
            % MNE fit
            params.modality = 'opm';
            params.chs = '*bz';
            fit_mne(save_path, opm_timelockedT, headmodel, sourcemodel, sourcemodel_inflated, params);
    
            if squid
                %% SQUID
                clear squid_timelocked        
                squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked;
                
                for i = 1:length(squid_timelocked)
                    squid_timelocked{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_squid.mat'])).squid_RS_cov;
                    if exist(fullfile(save_path, [params.sub '_ER_squid.mat']),'file')
                        squid_timelocked{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_squid.mat'])).squid_ER_cov;
                    end
                end
        
                params.modality = 'squidgrad';
                params.chs = 'meggrad';
                fit_mne(save_path, squid_timelocked, headmodel, sourcemodel, sourcemodel_inflated, params);
            end
        end
    end
end

save(fullfile(base_save_path, 'group_results.mat'), 'grp_tag_opm','grp_tag_squid','grp_SNR_opm','grp_SNR_squid','grp_pp_opm','grp_pp_squid');

%% clear and close all, then exit to free memory
close all
clear all
exit

%% Functions
function [subjects, sessions] = getSubjectsAndSessions(folderPath,natmeg)
    % Get a list of all subject folders
    if natmeg
        subjectFolders = dir(fullfile(folderPath, 'NatMEG_*'));
    else
        subjectFolders = dir(fullfile(folderPath));
    end
    % Initialize cell arrays for subjects and sessions
    subjects = {};
    sessions = {};
    
    i_sub = 0;
    for i = 1:length(subjectFolders)
        if subjectFolders(i).isdir && length(subjectFolders(i).name)>2
            i_sub = i_sub + 1;
            subjects{i_sub,1} = subjectFolders(i).name;

            sessionFolders = dir(fullfile(folderPath, subjectFolders(i).name));
            i_ses = 0;
            for j = 1:length(sessionFolders)
                if sessionFolders(j).isdir && length(sessionFolders(j).name) == 6 && all(isstrprop(sessionFolders(j).name, 'digit'))
                    i_ses = i_ses + 1;
                    sessions{i_sub, i_ses} = sessionFolders(j).name;
                end
            end
        end
    end
end
