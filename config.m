%% This config file contains functions to be called by the main.m scripts.
function out = config(what, varargin)
    % Parse optional input (e.g. modality = 'opm', 'meg', or 'eeg')
    switch what
        case 'params'
            if ~isempty(varargin)
                modality = varargin{1};
            else
                modality = 'meg';  % default modality
            end
            out = get_params(modality);
        case 'overwrite'
            out = get_overwrite();
        case 'paradigm'
            out = get_paradigm();
        case 'paths'
            out = get_paths();
        case 'skip'
            out = get_skip();
    end
end


%% Overwrite specific steps in the analysis 
function overwrite = get_overwrite()
    overwrite = [];
    overwrite.preproc = true;
    overwrite.timelock = true;
    overwrite.mri = false;
    overwrite.coreg = false;
    overwrite.dipole = false;
    overwrite.mne = false;
end

%% Params 
function params = get_params(modality)
    params = [];
    % Trial lenght
    params.pre = 0.3; % Trial prestim in seconds
    params.post = 0.8; % Trial poststim in seconds
    params.pad = 0.2; % Trial (pre and post) padding in seconds
    params.delay = 0.01; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).

    % Filter settings
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

    % ICA settings
    params.n_comp = 40; % Number of ICA components
    params.manual_ica = false; % Manually select ICA components to remove?
    params.save_ica = 1; % Save plots and components
    params.ica_cor = 0.8; % Cutoff for correlation with EOG/ECG 
    params.ica_coh = 0.95; % Cutoff for coherence with EOG/ECG 

    % BADS detection settings
    params.corr_threshold = 0.7; % Correlation threshold for badchannel neighbors
    params.z_threshold = 20; % Zmax threshold for badchannel and trial detection
    params.opm_std_threshold = 5e-12; % Stddev threshold for badtrial detection
    params.squid_std_threshold = 2.5e-12; % Stddev threshold for badtrial detection

    % HPI settings
    params.hpi_freq = 33; % HPI coil frequency
    params.hpi_gof = 0.9; % Minimum goodness-of-fit for including coil in hpi analysis

    % Source model settings
    params.src_density = '8'; % Sourcemodel density ('4', '8' or '32') = approximate number of sources per hemisphere
    params.source_fixedori = true; % use fixed orientation sources (along vertex normals); if false: use three orthogonal sources per location
    params.use_cov = 'resting_state'; % noise cov to use; default= ' ' for prestim, alt: 'resting_state', 'empty_room'

    % Peak settings
    params.peaks = []
    params.peaks.label = {'M100', 'M300'}; % Peak labels
    params.peaks.peak_latency = {0.1, 0.3}; % Peak latencies in seconds

    switch modality
        case 'opm'
            params.modality = 'opm';
            params.layout = 'fieldlinebeta2bz_helmet.mat'; % Layout for OPM data
            params.chs = '*_b*'; % Channel pattern for OPM data
            params.amp_scaler = 1e15;
            params.amp_label = 'B [fT]';
        case 'meg'
            params.modality = 'squid';
            params.layout = 'neuromag306mag.lay';
            params.chs = 'meg';
    end
end

function paradigm = get_paradigm()
    paradigm.paradigms = {'AudOdd'}; % Paradigms to analyze for all participants and sessions
    paradigm.trigger_codes = {3 5 11 13}; % Trigger values to timelock
    paradigm.trigger_labels = {'low_nogo', 'low_go', 'high_nogo', 'high_go'}; % Labels corresponding to the trigger values    
end

function paths = get_paths()
    paths = [];
    paths.base_data_path = '~/../../projects/capsi/bids';
    paths.base_save_path = '~/temp_output';
    paths.base_matlab_path = '/usr/local/MATLAB/R2024b';
    paths.project_scripts_path = '~/opm_general';
end

function skip = get_skip()
    skip = [];
    skip.subjects = {'sub-001', 'sub-002', 'sub-003'}; % Subjects to skip entirely (make sure they are zero padded to the amount of digitis in your bids foder)
    skip.sessions = {}; % Sessions to skip entirely
    skip.subsessions = {'sub-001_ses-01', 'sub-002_ses-02'}; % Subjects-session combinations to skip 
end