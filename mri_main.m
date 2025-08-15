%% Reset all
clear all
close all
restoredefaultpath

% make sure that pwd is opm_general in the users home folder
disp(pwd)
addpath(pwd)
overwrite = config('overwrite');
params = config('params', 'opm');
paradigm = config('paradigm');
paths = config('paths');
skip = config('skip');

%% Set up fieldtrip (now assumed to be in users home folder)
addpath(fullfile(pwd,'../fieldtrip')) % Fieldtrip path
ft_defaults

%% Subjects + dates
subsessions = readtable(fullfile(pwd, '../completed_meg_sessions_long_DATA_2025-07-23.csv'), 'Delimiter',',');
subject_list = unique(subsessions.new_subject_id);

mri_files = {'00000001.dcm' 
    '/mri/sub-15931_T1w.nii.gz'  
    '/nifti/anat/sub-15985_T1w.nii.gz'};

%% Loop over subjects
for i_sub = 1:length(subject_list)
    params.sub = ['sub-' num2str(subject_list(i_sub),'%03d')];

    sessions = subsessions(subsessions.new_subject_id==subject_list(i_sub),:).new_session_id;
    dates = subsessions(subsessions.new_subject_id==subject_list(i_sub),:).old_session_id;

    natmeg_id = unique(subsessions(subsessions.new_subject_id==subject_list(i_sub),:).old_subject_id);

    for i_ses = 1:length(sessions)
        % create zeropadded sessions
        params.ses = ['ses-' num2str(sessions(i_ses),'%02d')];
        %% Paths
        % the prepare_mri function needs a path to MRI and a path to the file containing the polhemus headshape (aux_file)
        %aux_path = fullfile(paths.base_data_path, params.sub, params.ses, 'eeg');

        % currently have to use eeg.fif from RAW instead of bids
        aux_path = fullfile('~/../../projects/capsi/raw/squid',['NatMEG_' num2str(natmeg_id)], num2str(dates(i_ses)),'meg');

        % similarly get MRI data from raw instead of bids
        mri_path = fullfile(paths.base_raw_path,'mri',params.sub, params.ses);
        save_path_mri = fullfile(paths.base_save_path, params.sub, params.ses, 'mri');
        
        % Create folders if they do not exist yet
        if ~exist(save_path_mri, 'dir')
            mkdir(save_path_mri)
        end

        % get correct aux file
        aux_files = [];
        tmp_eeg = dir(fullfile(aux_path, ['*' 'Aud' '*' 'EEG' '*' '.fif']));

        try
            aux_file = fullfile(tmp_eeg.folder, tmp_eeg.name);
        catch
            disp('No auxiliary files found');
            disp(params.sub);
            disp(params.ses);
            continue
        end

        %mri_file = fullfile(mri_path, 'orig','001.mgz');
        prepare_mri(mri_path,aux_file,save_path_mri,params);

        disp(mri_path)
        disp(aux_file)
        disp(save_path_mri)
        close all
    end
end