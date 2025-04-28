%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/23106_opmbci/';
    base_save_path = '/home/chrpfe/Documents/23106_opmbci';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/23106_opmbci/opmbci_preprocessing';
else
    % Laptop:
    base_data_path = '/Volumes/dataarchvie/23106_opmbci';
    base_save_path = '/Users/christophpfeiffer/data_local/23106_opmbci';
    base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/christophpfeiffer/opm_general';
end

%% Set up fieldtrip
addpath(fullfile(base_matlab_path,'fieldtrip-20231220/')) % Fieldtrip path
addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Overwrite
overwrite = [];
overwrite.read = true;
overwrite.preproc = true;
overwrite.mri = false;
overwrite.coreg = false;
overwrite.sourcerec = false;

%% Params
params = [];
params.pre = 0.2; %sec
params.post = 1.4; %sec
params.pad = 0.2;
params.filter = [];
params.filter.hp_freq = 1;
params.filter.lp_freq = 100;
params.filter.bp_freq = [];
params.filter.notch = sort([50:50:150 60:60:120]);
params.n_comp = 40;
params.manual_ica = false;
params.ica_cor = 0.8; % cutoff for EOG/ECG coherence
params.ica_coh = 0.95; % cutoff for EOG/ECG coherence
params.save_ica = 1; % save plots and components
params.z_threshold = 20;
params.corr_threshold = 0.7; % correlation threshold for badchannel neighbors
params.opm_std_threshold = 2.5e-12;
params.hpi_freq = 33;
params.hpi_gof = 0.9;

params.trigger_codes = [1 2 3 4 5 6];
params.trigger_labels = {'D', 'P', 'T', 'F', 'O', 'G'};

params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';

%% Subjects + dates
sub = {'NatMEG_0953'};

ses = {'250203'};

paradigms = {'RSEO'; 'RSEC'};

%%
i_sub = 1;
for i_ses = 1:length(ses)
    %% Loop over subjects
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    params.ses = ['ses_' num2str(i_ses,'%02d')];
    
    %% Paths
    raw_path = fullfile(base_data_path, sub{i_sub}, ses{i_ses});
    save_path = fullfile(base_save_path, params.sub, params.ses);
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end
    if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
    end
    for i_paradigm = 1:length(paradigms)
        opm_files = fullfile(raw_path,'osmeg',[paradigms{i_paradigm} 'OPM_raw.fif']); % opm files 
        aux_files = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'EEG.fif']); % corresponding aux files containing EOG/ECG
    end
    hpi_path = fullfile(raw_path,'osmeg');
    mri_path = '/Volumes/dataarchvie/23106_opmbci/NatMEG_0953/mri/';
    mri_file = fullfile(mri_path,'mri','orig','001.mgz');
    
    for i_paradigm = 1:length(paradigms)
        params.paradigm = paradigms{i_paradigm};
        
        %% Read and preproc
        if overwrite.read == true
            % Read data and reject bad channels
            disp(['Reading file: ' num2str(i_file) '/' num2str(length(i_paradigm)) '...'])
            ft_hastoolbox('mne', 1);
            data = read_osMEG(opm_files{i_file}, aux_files{i_file}, save_path, params); % Read data
            save(fullfile(save_path, [params.sub '_' params.paradigm '_data']), 'data',"-v7.3"); disp('done');
        end

        if overwrite.preproc == true
            % ICA
            disp('Running ICA ...')
            data_ica = ica_MEG(data, save_path, params);
            save(fullfile(save_path, [params.sub '_' params.paradigm '_data_ica']), 'data_ica',"-v7.3"); disp('done');
        end
    end
    
    %% MRI
    ft_hastoolbox('mne',1);
    if overwrite.mri==true
        ft_hastoolbox('mne', 1);
        prepare_mri_opmbci(mri_file,aux_files(1),save_path)

        % Read and transform cortical restrained source model
        files = dir(fullfile(mri_path,'workbench'));
        for i = 1:length(files)
            if endsWith(files(i).name,'.L.midthickness.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            end
            if endsWith(files(i).name,'.L.aparc.8k_fs_LR.label.gii')
                filename2 = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        clear sourcemodel mri_resliced
        mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
        sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
        
        aparc_L = ft_read_atlas({filename2,filename});
        aparc_R = ft_read_atlas({strrep(filename2,'.L.','.R.'),strrep(filename,'.L.','.R.')});
        tmp = ft_read_atlas(strrep(filename2, '.L.', '.R.'),'format','caret_label');
        n_labels = length(aparc_L.parcellationlabel);
        atlas = [];
        atlas.parcellationlabel = [aparc_L.parcellationlabel; aparc_R.parcellationlabel];
        atlas.parcellation = [aparc_L.parcellation; aparc_R.parcellation + n_labels];
        atlas.rgba = [aparc_L.rgba; aparc_R.rgba; [0 0 0 1]];
        n_labels = length(atlas.parcellationlabel);
        atlas.parcellation(isnan(atlas.parcellation))=n_labels+1;
        sourcemodel.brainstructure = atlas.parcellation;
        sourcemodel.brainstructurelabel = atlas.parcellationlabel;
        sourcemodel.brainstructurecolor = atlas.rgba;
    
        T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
        sourcemodel = ft_transform_geometry(T, sourcemodel);
        sourcemodel.inside = true(size(sourcemodel.pos,1),1);
        save(fullfile(save_path, [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
    end
    
    %% HPI localization
    ft_hastoolbox('mne',1);
    
    if overwrite.coreg==true
        ft_hastoolbox('mne', 1);
        params.include_chs = load(fullfile(save_path, ['include_chs' num2str(length(opm_files))])).include_chs;
        %data_ica = load(fullfile(save_path, [params.sub '_' params.modality '_auditory' num2str(i_file) '.mat'])).data_ica;
        %params.include_chs = data_ica.label;   
        clear data_ica
        opm_trans = fit_hpi_opmbci(hpi_path, aux_files{1}, save_path, params);
    
        for i_file = 1:length(opm_files)
            clear timelocked
            data_ica = load(fullfile(save_path, [params.sub '_' params.modality '_auditory' num2str(i_file) '.mat'])).data_ica;
            data_ica.grad.chanpos = opm_trans.transformPointsForward(data_ica.grad.chanpos);
            data_ica.grad.coilpos = opm_trans.transformPointsForward(data_ica.grad.coilpos);
            data_ica.grad.chanori = (opm_trans.Rotation'*data_ica.grad.chanori')';
            data_ica.grad.coilori = (opm_trans.Rotation'*data_ica.grad.coilori')';
            save(fullfile(save_path, [params.sub '_' params.modality '_auditoryT' num2str(i_file)]),'data_ica',"-v7.3");
            cfg = [];
            cfg.covariance = 'yes';
            cfg.covariancewindow = [-params.pre 0];
            timelocked = ft_timelockanalysis(cfg,data_ica);
            save(fullfile(save_path, [params.sub '_' params.modality '_auditory_timelockedT' num2str(i_file)]),'timelocked',"-v7.3");
            clear data_ica timelocked
        end

        % Plot source and head models
        clear headmodels sourcemodel
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
        sourcemodel = load(fullfile(save_path, 'sourcemodel.mat')).sourcemodel;
    
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(timelocked.grad,'unit','cm')
        hold off;
        title('OPM-MEG')
        view([-140 10])
        saveas(h, fullfile(save_path, 'figs', 'opm_layout.jpg'))
        close all
    
        clear timelocked mri_resliced opm_trans
    end
    
    %% MNE
    ft_hastoolbox('mne',1);
    if overwrite.sourcerec==true
        clear headmodels sourcemodel
        sourcemodel = load(fullfile(save_path, 'sourcemodel')).sourcemodel;
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
        pos = sourcemodel.pos;
        tri = sourcemodel.tri;
        roi = sourcemodel.brainstructurelabel;
        
        N_rois = length(sourcemodel.brainstructurelabel);
        N_sources = size(pos, 1);
        mapping_matrix = zeros(N_rois, N_sources);
        for i = 1:N_rois
            mapping_matrix(i,sourcemodel.brainstructure==i) = 1;
        end
        roi_counts = sum(mapping_matrix, 2);
        mapping_matrix = mapping_matrix ./ repmat(roi_counts,[1 size(mapping_matrix,2)]);
    
        for i_file = 1:length(opm_files)
            if i_ses == 1
                data_set = '';
            else
                data_set = num2str(i_file);
            end
        
            timelocked = load(fullfile(save_path, [params.sub '_' params.modality '_motorimag_timelockedT' data_set '.mat'])).timelocked;
            mne = fit_mne_opmbci(timelocked,headmodels,sourcemodel,params); 
            mne_inv = zeros(length(mne.avg.filter),length(mne.avg.label));
            for i = 1:length(mne.avg.filter)
                mne_inv(i,:) = mne.avg.filter{i};
            end
            % Plot timelocked
            h = figure ;
            plot(timelocked.time,timelocked.avg)
            title('Motor imagery timelocked')
            saveas(h, fullfile(save_path, 'figs', 'opm_motorimag_timelocked.jpg'))
            close all
            clear mne timelocked
            data_ica = load(fullfile(save_path, [params.sub '_' params.modality '_motorimag' data_set '.mat'])).data_ica;
            trial = cell(length(data_ica.trial),1);
            for i_trl = 1:length(data_ica.trial)
                trial{i_trl} = mapping_matrix*((mne_inv*data_ica.trial{i_trl}));
            end
            time = data_ica.time{1};
            label = params.trigger_labels(data_ica.trialinfo-64);
            save(fullfile(save_path, ['motorimag_mne' num2str(i_file) '.mat']),'mne_inv','trial','time','label','pos','tri','roi','-v7.3');
            clear mne_inv
        end
    end
    close all
end

%% clear and close all, then exit to free memory
close all
clear all
exit

%% Functions
function files = findOpmFiles(directory, pattern)
    % This function finds all files in the specified directory that match
    % the pattern "TrainingSet" followed by a number from 1 to 10 and then "_raw.fif".
    
    % Define the pattern
    pattern = '*_raw.fif';
    
    % Get a list of all files in the directory
    allFiles = dir(directory);
    
    % Initialize an empty cell array to store matching files
    files = {};
    
    % Loop through all files and check if they match the pattern
    for i = 1:length(allFiles)
        if ~allFiles(i).isdir
            if ~isempty(regexp(allFiles(i).name, pattern, 'once'))
                files{end+1} = fullfile(directory, allFiles(i).name); %#ok<AGROW>
            end
        end
    end
end