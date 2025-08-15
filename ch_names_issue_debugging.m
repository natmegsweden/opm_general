%% Set up
clear all
close all
restoredefaultpath

% make sure that pwd is opm_general in the users home folder
disp(pwd)
addpath(pwd)

% Get params from config
overwrite = config('overwrite');
params = config('params', 'opm');
paradigm = config('paradigm');
paths = config('paths');
skip = config('skip');

% Set up fieldtrip (now assumed to be in users home folder)
addpath(fullfile(pwd,'../fieldtrip')) % Fieldtrip path
ft_defaults

% Subjects + dates
subsessions = readtable(fullfile(pwd, '../completed_meg_sessions_long_DATA_2025-07-23.csv'), 'Delimiter',',');
subject_list = unique(subsessions.new_subject_id);


%% Pick a subject
i_sub = 3;
i_ses = 2;
i_paradigm = 1;

params.sub = ['sub-' num2str(subject_list(i_sub),'%03d')];
disp(params.sub)


%%
save_path = fullfile(paths.base_save_path, params.sub);  
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

% find sessions for subject
sessions = subsessions(subsessions.new_subject_id==subject_list(i_sub),:).new_session_id;
dates = subsessions(subsessions.new_subject_id==subject_list(i_sub),:).old_session_id;
natmeg_id = unique(subsessions(subsessions.new_subject_id==subject_list(i_sub),:).old_subject_id);

disp(natmeg_id)

if length(natmeg_id)>1
    error('Multiple natmeg ids for subject %s', params.sub)
end

params.ses = ['ses-' num2str(sessions(i_ses),'%02d')];

fprintf('working on %s %s \n',params.sub,params.ses);

%% Get subject-session paths
raw_path = fullfile(paths.base_data_path, params.sub, params.ses);
save_path = fullfile(paths.base_save_path, params.sub, params.ses);

% temporarily get auxilary data from raw instead of bids folder and also read OPM data from raw folder
tmp_eeg_raw_path = fullfile('~/../../projects/capsi/raw/squid',['NatMEG_' num2str(natmeg_id)], num2str(dates(i_ses)),'meg');
tmp_opm_raw_path = fullfile('~/../../projects/capsi/raw/opm',['sub-' num2str(natmeg_id)]);
opm_files = [];
% search for AudOdd data, might be split
tmp = dir(fullfile(raw_path, 'meg',['*' paradigm.paradigms{i_paradigm} '_acq-hedscan' '*' '_meg.fif']));

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
tmp_eeg = dir(fullfile(tmp_eeg_raw_path, [paradigm.paradigms{i_paradigm} '*' 'EEG' '*' '.fif']));
tmp_opm = dir(fullfile(tmp_opm_raw_path, ['*' paradigm.paradigms{i_paradigm} '*' '.fif']));
aux_files{i_paradigm} = fullfile(tmp_eeg.folder, tmp_eeg.name); % corresponding aux files containing EOG/ECG

if numel(tmp_opm) > 1
    % Return a cell array of full paths
    raw_opm_files{i_paradigm} = arrayfun(@(f) fullfile(f.folder, f.name), tmp_opm, 'UniformOutput', false);
else
    % Return a single path
    raw_opm_files{i_paradigm} = fullfile(tmp_opm.folder, tmp_opm.name);
end

% select session (date was not correct so no earlier filtering)
raw_opm_files{i_paradigm} = raw_opm_files{i_paradigm}(i_ses);

%% read BIDS data
cfg = [];
cfg.dataset        = opm_files{i_paradigm};
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_bids = ft_preprocessing(cfg);

%% read RAW data

cfg = [];
cfg.dataset        = raw_opm_files{i_paradigm};
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_raw = ft_preprocessing(cfg);

%% Channel names BIDS
bids_channel_names = opm_bids.label;
bid_orig_channel_names = opm_bids.hdr.orig{1}.orig.ch_names';

disp(bids_channel_names);
disp(bid_orig_channel_names);

%% Channel names raw
raw_channel_names = opm_raw.label;
raw_orig_channel_names = opm_raw.hdr.orig{1}.orig.ch_names';

disp(raw_channel_names);
disp(raw_orig_channel_names);
