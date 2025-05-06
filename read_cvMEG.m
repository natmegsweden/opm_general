function data = read_cvMEG(squid_file, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

if ~exist(squid_file,'file')
    error(['Did not find data file: ' squid_file])
end

%% --- Read triggers ---
trl = [];
cfg             = [];
cfg.datafile    = squid_file;
squid_raw         = ft_preprocessing(cfg);
squid_trig = find(contains(squid_raw.label,'STI101'));
trig = squid_raw.trial{1}(squid_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl(:,1) = find(trig)-(params.pre+params.pad)*squid_raw.fsample;
trl(:,2) = find(trig)+(params.post+params.pad)*squid_raw.fsample;
trl(:,3) = -(params.pre+params.pad)*squid_raw.fsample;
trl(:,4) = squid_raw.trial{1}(squid_trig,trig);
trl(:,1:2) = trl(:,1:2) + floor(0.041*squid_raw.fsample); % adjust for stim delay
trl = round(trl);

%% Data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
data = ft_preprocessing(cfg, data_raw);

cfg = [];
cfg.trl             = trl;
data = ft_redefinetrial(cfg,data);

cfg = [];
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
data = ft_preprocessing(cfg, data);

%% Remove padding
cfg = [];
cfg.latency = [-params.pre params.post];
data = ft_selectdata(cfg, data); 

%% Convert to sensor definitions to cm
data.grad = ft_convert_units(data.grad,'cm');

clear data_raw

end