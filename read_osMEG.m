function [data, badchs] = read_osMEG(opm_file, aux_file, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

if isempty(aux_file)
    opm_only = true;
else
    opm_only = false;
end

%% --- Read triggers ---
% OPM
trl_opm=[];
cfg = [];
cfg.datafile        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_raw = ft_preprocessing(cfg);
opm_trig = find(contains(opm_raw.label,'di'));
trig = opm_raw.trial{1}(opm_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_opm(:,1) = find(trig)-(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,2) = find(trig)+(params.post+params.pad)*opm_raw.fsample;
trl_opm(:,3) = -(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,4) = opm_raw.trial{1}(opm_trig,trig);
trl_opm(:,1:2) = trl_opm(:,1:2) + floor(0.041*opm_raw.fsample); % adjust for stim delay
trl_opm = round(trl_opm);

if ~opm_only
    % AUX
    trl_aux=[];
    cfg = [];
    cfg.datafile        = aux_file;
    aux_raw = ft_preprocessing(cfg);
    aux_trig = find(contains(aux_raw.label,'STI101'));
    trig = aux_raw.trial{1}(aux_trig,:)>0.5;
    trig = [false trig(2:end)&~trig(1:end-1)];
    trl_aux(:,1) = find(trig)-(params.pre+params.pad)*aux_raw.fsample;
    trl_aux(:,2) = find(trig)+(params.post+params.pad)*aux_raw.fsample;
    trl_aux(:,3) = -(params.pre+params.pad)*aux_raw.fsample;
    trl_aux(:,4) = aux_raw.trial{1}(aux_trig,trig);
    trl_aux(:,1:2) = trl_aux(:,1:2) + floor(0.041*aux_raw.fsample); % adjust for stim delay
    trl_aux = round(trl_aux);
    
    % Check if uneven amount of trial. If so assume error in beginning.
    if size(trl_aux,1) > size(trl_opm,1)
        trl_aux = trl_aux((end-size(trl_opm,1)+1):end,:);
    elseif size(trl_aux,1) < size(trl_opm,1)
        trl_opm = trl_opm((end-size(trl_aux,1)+1):end,:);
    end
    if trl_aux(:,4) ~= trl_opm(:,4) % Throw error if trials don't match.
        error('events do not match')
    end
    
    %% AUX data filter & epoch
    cfg = [];
    %cfg.datafile        = aux_file;
    %cfg.trl             = trl_aux;
    cfg.lpfilter        = 'yes';         
    cfg.lpfreq          = params.filter.lp_freq;
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
    %cfg.padding         = params.pre + params.post + 1;
    %cfg.paddingtype     = 'data';
    aux_epo = ft_preprocessing(cfg,aux_raw);
    
    cfg = [];
    cfg.trl             = trl_aux;
    aux_epo = ft_redefinetrial(cfg,aux_epo);
    
    cfg = [];
    cfg.dftfilter       = 'yes';        
    cfg.dftfreq         = params.filter.notch;
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-params.pre 0];
    aux_epo = ft_preprocessing(cfg,aux_epo);
end

%% OPM data filter & epoch
cfg = [];
%cfg.datafile        = opm_file;
%cfg.coordsys        = 'dewar';
%cfg.coilaccuracy    = 0;
%cfg.trl             = trl_opm;
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
%cfg.padding         = params.pre + params.post + 3;
%cfg.paddingtype     = 'data';
opm_epo = ft_preprocessing(cfg, opm_raw);

cfg = [];
cfg.trl             = trl_opm;
opm_epo = ft_redefinetrial(cfg,opm_epo);

cfg = [];
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
opm_epo = ft_preprocessing(cfg, opm_epo);

% Resample
if ~opm_only
    cfg            = [];
    cfg.time = aux_epo.time;
    cfg.detrend    = 'no';
    opm_epo_ds = ft_resampledata(cfg, opm_epo);
else
    cfg            = [];
    cfg.resamplefs = 1000;
    cfg.detrend    = 'no';
    opm_epo_ds = ft_resampledata(cfg, opm_epo);
end

%% Combine data
if ~opm_only
    EOG_channels = find(contains(aux_epo.label,'EOG'));
    ECG_channels = find(contains(aux_epo.label,'ECG'));
    EEG_channels = find(contains(aux_epo.label,'EEG'));
    MISC_channels = find(contains(aux_epo.label,'MISC'));
    TRIG_channels = find(contains(aux_epo.label,'STI101'));
    include_channels = [EOG_channels; ECG_channels; EEG_channels; MISC_channels; TRIG_channels];
    
    data = opm_epo_ds; 
    data.elec = aux_epo.elec;
    data.time = aux_epo.time;
    data.label = [data.label; aux_epo.label(include_channels)];
    data.hdr.label = data.label;
    data.hdr.nChans = data.hdr.nChans + length(include_channels);
    data.hdr.chantype = [data.hdr.chantype; aux_epo.hdr.chantype(include_channels)];
    data.hdr.chanunit = [data.hdr.chanunit; aux_epo.hdr.chanunit(include_channels)];
    data.sampleinfo = aux_epo.sampleinfo;
    data.trialinfo = aux_epo.trialinfo;
    n_smpl = size(aux_epo.trial{1},2);
    for i = 1:length(data.trial)
        data.trial{i} = [data.trial{i}(:,1:n_smpl); aux_epo.trial{i}(include_channels,:)]; 
    end
else
    data = opm_epo_ds;
end

%% Find bad opm channels
[badchs, badchs_flat, badchs_std, badchs_neighbors, badchs_outlier] = opm_badchannels(opm_raw, trl_opm, params);
save(fullfile(save_path, [params.paradigm '_badchs']), ...
    'badchs_flat', ...
    'badchs_std', ...
    'badchs_neighbors', ...
    'badchs_outlier',"-v7.3"); 
clear opm_raw

end