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
data_raw         = ft_preprocessing(cfg);
i_trig = find(contains(data_raw.label,'STI101'));
trig = data_raw.trial{1}(i_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl(:,1) = find(trig)-(params.pre+params.pad)*data_raw.fsample;
trl(:,2) = find(trig)+(params.post+params.pad)*data_raw.fsample;
trl(:,3) = -(params.pre+params.pad)*data_raw.fsample;
trl(:,4) = data_raw.trial{1}(i_trig,trig);
trl(:,1:2) = trl(:,1:2) + floor(params.delay*data_raw.fsample); % adjust for stim delay
trl = round(trl);

%% Data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
if ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
    if params.filter.hp_freq<1
        cfg.hpfilttype = 'firws';
    end
end
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


%% Reject jump trials
cfg = [];
cfg.channel = {'megmag'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_jump] = ft_badsegment(cfg, data);
data = ft_rejectartifact(cfg,data);

%% Reject noisy trials
cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'std';
cfg.threshold = params.squidmag_std_threshold;
[cfg,badtrl_squidmag_std] = ft_badsegment(cfg, data);
data = ft_rejectartifact(cfg,data);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'std';
cfg.threshold = params.squidgrad_std_threshold;
[cfg,badtrl_squidgrad_std] = ft_badsegment(cfg, data);
data = ft_rejectartifact(cfg,data);

%% Remove bad trials
[~,idx]=ismember(data.sampleinfo,badtrl_jump,'rows');
badtrl_jump = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidmag_std,'rows');
badtrl_std = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidgrad_std,'rows');
badtrl_std = unique([badtrl_std; find(idx)]);
save(fullfile(save_path, [params.paradigm '_badtrls']), ...
    'badtrl_jump', ...
    'badtrl_std', "-v7.3"); 

%% Remove padding
cfg = [];
cfg.latency = [-params.pre params.post];
data = ft_selectdata(cfg, data); 

%% Convert to sensor definitions to cm
data.grad = ft_convert_units(data.grad,'cm');

clear data_raw

end