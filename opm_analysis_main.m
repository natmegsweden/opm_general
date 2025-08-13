%% This main script runs the preprecessing of OPM data
% TO DO:
% - Figure out where completed_meg_sessions_long_DATA are stored
% - Figure out outputpath (now saved to ~/temp_output)
% - After proper BIDSification, use EEG.fif from bids structure, 
%       now uses from raw as the brainvision format does not have the head model 

%% Reset all
clear all
close all
restoredefaultpath

% make sure that pwd is opm_general in the users home folder
disp(pwd)
addpath(pwd)

%% Get params from config
overwrite = config('overwrite');
params = config('params', 'opm');
paradigm = config('paradigm');
paths = config('paths');
skip = config('skip');

%% Set up fieldtrip (now assumed to be in users home folder)
addpath(fullfile(pwd,'../fieldtrip')) % Fieldtrip path
ft_defaults

%% Analyse squid data?
squid = false;

%% Subjects + dates
subsessions = readtable(fullfile(pwd, '../completed_meg_sessions_long_DATA_2025-07-23.csv'), 'Delimiter',',');
subject_list = unique(subsessions.new_subject_id);
natmeg_ids = unique(subsessions.old_subject_id);

if ~length(subject_list)==length(natmeg_ids)
    error('Mismatch between new and old subject IDs')
end

%% Loop over subjects
for i_sub = 1:length(subject_list)
    % zeropad to three zeros
    params.sub = ['sub-' num2str(subject_list(i_sub),'%03d')];
    disp(params.sub)
    disp(natmeg_ids(i_sub))

    % check if subject is in skip list
    if ismember(params.sub, skip.subjects)
        disp(['Skipping subject: ' params.sub])
        continue
    end
    
    % create save path
    save_path = fullfile(paths.base_save_path, params.sub);  
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end

    % find sessions for subject
    sessions = subsessions(subsessions.new_subject_id==subsessions.new_subject_id(i_sub),:).new_session_id;
    dates = subsessions(subsessions.new_subject_id==subsessions.new_subject_id(i_sub),:).old_session_id;
    
    %% Loop over sessions
    for i_ses = 1:length(sessions)
        % create zeropadded sessions
        params.ses = ['ses-' num2str(sessions(i_ses),'%02d')];

        % skip if session in skip.sessions or combination of subject and session in skip.subsessions
        if ismember(params.ses, skip.sessions) || ismember([params.sub '_' params.ses], skip.subsessions)
            disp(['Skipping session: ' params.ses ' for subject: ' params.sub])
            continue
        end

        fprintf('working on %s %s \n',params.sub,params.ses);

        %% Get subject-session paths
        raw_path = fullfile(paths.base_data_path, params.sub, params.ses);
        save_path = fullfile(paths.base_save_path, params.sub, params.ses);

        % temporarily get auxilary data from raw instead of bids folder
        tmp_eeg_raw_path = fullfile('~/../../projects/capsi/raw/squid',['NatMEG_' num2str(natmeg_ids(i_sub))], num2str(dates(i_ses)),'meg');
        if ~exist(save_path, 'dir')
           mkdir(save_path)
        end
        if ~exist(fullfile(save_path,'figs'), 'dir')
           mkdir(fullfile(save_path,'figs'))
        end
        for i_paradigm = 1:length(paradigm.paradigms)
            opm_files = [];
            % search for AudOdd data, might be split
            tmp = dir(fullfile(raw_path, 'meg',['*' paradigm.paradigms{i_paradigm} '_acq-hedscan' '*' '_meg.fif']));
            if isempty(tmp)
                disp(['No OPM data found for paradigm: ' paradigm.paradigms{i_paradigm} ' in ' fullfile(raw_path, 'meg')])
                continue
            else
                %disp(tmp)
                % handle split files
                 if numel(tmp) > 1
                    % Return a cell array of full paths
                    opm_files{i_paradigm} = arrayfun(@(f) fullfile(f.folder, f.name), tmp, 'UniformOutput', false);
                else
                    % Return a single path
                    opm_files{i_paradigm} = fullfile(tmp.folder, tmp.name);
                end
                disp(opm_files{i_paradigm})

                % temporarily look for EEG.fif in raw folder instead of brainvision eeg.eeg in bids folder
                %tmp_eeg = dir(fullfile(raw_path,'eeg',['*' paradigm.paradigms{i_paradigm} '_acq-triux' '*' '_eeg.eeg']));
                tmp_eeg = dir(fullfile(tmp_eeg_raw_path, [paradigm.paradigms{i_paradigm} '*' 'EEG.fif']));
                aux_files{i_paradigm} = fullfile(tmp_eeg.folder, tmp_eeg.name); % corresponding aux files containing EOG/ECG
            end
        end
        if isempty(opm_files)
            disp(['No OPM data found for subject: ' params.sub ' in session: ' params.ses])
            continue
        end
        hpi_path = fullfile(raw_path,'meg');
        %% Loop over paradigm.paradigms/tasks
        for i_paradigm = 1:length(paradigm.paradigms)
            params.paradigm = paradigm.paradigms{i_paradigm};
            params.trigger_codes = paradigm.trigger_codes; 
            params.trigger_labels = paradigm.trigger_labels; 
            
            disp(['Processing paradigm: ' params.paradigm])
            %% Read and preproc
            if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica.mat']),'file')
                ft_hastoolbox('mne', 1);
    
                % Read data 
                disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigm.paradigms)) '...'])
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
                timelocked = timelock(data_ica, save_path, params);
                save(fullfile(save_path, [params.paradigm '_timelocked']), 'timelocked', '-v7.3'); 
                clear timelocked
            end
            clear data_ica
    
            % if squid
            %     %% Read and preproc - SQUID-MAG
            %     params.modality = 'squid';
            %     params.layout = 'neuromag306mag.lay';
            %     params.chs = 'meg';
        
            %     if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat']),'file')
            %         ft_hastoolbox('mne', 1);
        
            %         % Read data 
            %         disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigm.paradigms)) '...'])
            %         data_epo = read_cvMEG(squid_files{i_paradigm}, params); % Read data
                    
            %         % ICA
            %         disp('Running ICA ...')
            %         if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1 % No ExG data
            %             params.manual_ica = 1;
            %             params.save_ica = 1;
            %         end
            %         data_ica = ica_MEG(data_epo, save_path, params);
            %         save(fullfile(save_path, [params.paradigm '_data_ica_squidmag']), 'data_ica',"-v7.3"); disp('done');
            %         clear data_epo
            %     else
            %         data_ica = load(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat'])).data_ica;
            %     end
                
            %     if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
            %         params.modality = 'squidmag';
            %         params.layout = 'neuromag306mag.lay';
            %         params.chs = 'megmag';
            %         params.amp_scaler = 1e15;
            %         params.amp_label = 'B [fT]';
            %         timelocked = timelock(data_ica, save_path, params);
            %         save(fullfile(save_path, [params.paradigm '_timelocked_squidmag']), 'timelocked', '-v7.3'); 
            %         clear timelocked
            %     end
            %     clear data_ica
            % end
        end
        
        %% HPI localization
        ft_hastoolbox('mne',1);
        
        % create save path for mri and check if folder exists, if not, create
        save_path_mri = fullfile(save_path, 'mri');
        if ~exist(save_path_mri, 'dir')
            mkdir(save_path_mri)
        end
        
        if exist(fullfile(save_path_mri, 'opm_trans.mat'),'file') && overwrite.coreg==false
            disp(['Not overwriting OPM transform for ' params.sub]);
        else
            ft_hastoolbox('mne', 1);
            params.include_chs = load(fullfile(save_path, [params.paradigm '_include_chs.mat'])).include_chs;
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

save(fullfile(paths.base_save_path, 'group_results.mat'), 'grp_tag_opm','grp_tag_squid','grp_SNR_opm','grp_SNR_squid','grp_pp_opm','grp_pp_squid');

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
