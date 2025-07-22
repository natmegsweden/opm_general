%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/21099_opm/';
    base_save_path = '/home/chrpfe/Documents/21099_opm/';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/21099_opm/phalanges';
    on_server = true;
else
    % Laptop:
    base_data_path = '/Volumes/dataarchvie/21099_opm';
    base_save_path = '/Users/christophpfeiffer/data_local/Benchmarking/';
    base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/christophpfeiffer/opm_phalanges';
    on_server = false;
end

%% Set up fieldtrip
addpath(fullfile(base_matlab_path,'fieldtrip/')) % Fieldtrip path
addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Subjects + dates
subses = {'0005' '240208';
    '0905' '240229';
    '0916' '240320';
    '0953' '241104';
    '1096' '241022';
    '1153' '240321';
    '1167' '240425';
    '1186' '240925';
    '1190' '241023';
    '1191' '241024';
    '1193' '241029';
    '1194' '241029';
    '1195' '241030'};
mri_files = {'00000001.dcm' 
    '/mri/sub-15931_T1w.nii.gz'  
    '/nifti/anat/sub-15985_T1w.nii.gz'};

if on_server
    subs_to_run = 1:size(subses,1);
else
    subs_to_run = 2; %1:size(subses,1)
end
excl_subs = [1];

%% Loop over subjects
for i_sub = setdiff(subs_to_run,excl_subs)
    params.sub = ['sub_' num2str(i_sub,'%02d')];

    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    save_path_mri = fullfile(base_save_path,'MRI',params.sub);
    
    % Create folders if they do not exist yet
    if ~exist(fullfile(base_save_path,'MRI'), 'dir')
        mkdir(fullfile(base_save_path,'MRI'))
    end
    if ~exist(save_path_mri, 'dir')
        mkdir(save_path_mri)
    end
    meg_file = fullfile(raw_path, 'meg', 'AudOddMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
    %mri_file = fullfile(mri_path, 'orig','001.mgz');
    prepare_mri(mri_path,meg_file,save_path_mri);
    close all
end