function [data] = read_osMEG(opm_file, aux_file, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% opm_file: path to the file containing the OPM-MEG data
% aux_file: path to synchronously recorded file containing ECG, EOG and/or
%           EEG. If not appliccable set to []
% save_path: path where to save results
% params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

% since OPM data is sometimes split into multiple files, this function
% can handle opm_file being a cell array of file paths or a single file path. 

num_files = numel(opm_file);

if num_files > 1
    for i = 1:num_files
        if ~exist(opm_file{i},'file')
            error(['Did not find OPM file: ' opm_file{i}])
        end
    end
else
    if ~exist(opm_file,'file')
        error(['Did not find OPM file: ' opm_file])
    end
end

if isempty(aux_file)
    opm_only = true;
else
    opm_only = false;
    if ~exist(aux_file,'file')
         error(['Did not find AUX file: ' aux_file])
    end
end

%% --- Read triggers ---
% OPM
trl_opm=[];
cfg = [];
cfg.dataset        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_raw = ft_preprocessing(cfg);

% Correct channel names from header and remove -sXX sensor names to match grad structure 
opm_raw.label = opm_raw.hdr.orig{1}.orig.ch_names';
opm_raw.label = regexprep(opm_raw.label, '-s\d+$', '');

% Find trigger channel 
i_trig_opm = find(contains(opm_raw.label,'di'));

% hacky solution for now:

%{
 for i = 1:length(opm_raw.label)
    % check for channel that is 0 over first 10 samples but not the first 20000 samples to avoid picking HPI channels 
    if sum(opm_raw.trial{1}(i,1:10)) == 0 & sum(opm_raw.trial{1}(i,1:200000))
        disp(i)
        % set trigger opm to value
        i_trig_opm = i
        % Create plot
        fig = figure('Visible', 'off');  % Don’t show plot (faster)
        plot(opm_raw.trial{1}(i,:));
        title(sprintf('Sensor %d', i));

        % Create filename with padded index (e.g., plot_001.png)
        filename = sprintf('plot_%03d.png', i);

        % Save the figure
        saveas(fig, filename);

        % Close the figure to save memory
        close(fig);
    end
end 
%}

trig = opm_raw.trial{1}(i_trig_opm,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_opm(:,1) = find(trig)-(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,2) = find(trig)+(params.post+params.pad)*opm_raw.fsample;
trl_opm(:,3) = -(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,4) = opm_raw.trial{1}(i_trig_opm,trig);
trl_opm(:,1:2) = trl_opm(:,1:2) + floor(params.delay*opm_raw.fsample); % adjust for stim delay
trl_opm = round(trl_opm);

disp(['-----------i_trig_opm-----------'])
disp(i_trig_opm)
disp(trl_opm)
disp(['----------------------'])

if ~opm_only
    % AUX
    trl_aux=[];
    cfg = [];
    cfg.datafile        = aux_file;
    aux_raw = ft_preprocessing(cfg);
    i_trig_aux = find(contains(aux_raw.label,'STI101'));
    trig = aux_raw.trial{1}(i_trig_aux,:)>0.5;
    trig = [false trig(2:end)&~trig(1:end-1)];
    trl_aux(:,1) = find(trig)-(params.pre+params.pad)*aux_raw.fsample;
    trl_aux(:,2) = find(trig)+(params.post+params.pad)*aux_raw.fsample;
    trl_aux(:,3) = -(params.pre+params.pad)*aux_raw.fsample;
    % trl_aux(:,4) = aux_raw.trial{1}(i_trig_aux,trig)/1000000; % trigger values are 1x10^6 instead of 1 so add adjustment here when using brainvision format
    trl_aux(:,4) = aux_raw.trial{1}(i_trig_aux,trig); 
    trl_aux(:,1:2) = trl_aux(:,1:2) + floor(params.delay*aux_raw.fsample); % adjust for stim delay
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
    cfg.lpfilter        = 'yes';         
    cfg.lpfreq          = params.filter.lp_freq;
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
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
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
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
    %data.elec = aux_epo.elec; % elec is not in aux_epo
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

%% Find & remove bad opm channels
[badchs, badchs_flat, badchs_std, badchs_neighbors, badchs_outlier] = opm_badchannels(opm_raw, trl_opm, params);
save(fullfile(save_path, [params.paradigm '_badchs']), ...
    'badchs_flat', ...
    'badchs_std', ...
    'badchs_neighbors', ...
    'badchs_outlier',"-v7.3"); 

cfg = [];

% save good channels
include_chs = setdiff(data.label,badchs);
cfg.channel = include_chs;
save(fullfile(save_path, [params.paradigm '_include_chs']), ...
    'include_chs',"-v7.3");
data = ft_selectdata(cfg, data);
hdr = data.hdr;

%% Spatiotemporal filtering
cfg = []; % separate ExG channels
cfg.channel = {'EOG*', 'ECG*'};
ExG = ft_selectdata(cfg,data);

if params.apply_hfc
    cfg = [];
    cfg.channel = '*bz';
    cfg.order = params.hfc_order;
    cfg.residualcheck = 'no';
    data = ft_denoise_hfc(cfg, data);
elseif params.apply_amm
    cfg = [];
    cfg.channel = '*bz';
    cfg.updatesens = 'yes';
    cfg.residualcheck = 'no';
    cfg.amm = [];
    cfg.amm.order_in = params.amm_in;
    cfg.amm.order_out = params.amm_out;
    cfg.amm.thr = params.amm_thr;
    data = ft_denoise_amm(cfg, data);
else
    data = data;
end

% Recombine with ExG channels
data.label = vertcat(data.label,ExG.label);

data.hdr = hdr;
incl = ismember(hdr.label,data.label);
data.hdr.label = hdr.label(incl);
data.hdr.chantype = hdr.chantype(incl);
data.hdr.chanunit = hdr.chanunit (incl);
for i = 1:length(data.trial)
    data.trial{i} = vertcat(data.trial{i}, ExG.trial{i}); 
end

%% Reject jump trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_jump] = ft_badsegment(cfg, data);
data = ft_rejectartifact(cfg,data);

%% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'std';
cfg.threshold = params.opm_std_threshold;
[cfg,badtrl_std] = ft_badsegment(cfg, data);
data = ft_rejectartifact(cfg,data);

%% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) && params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    data = ft_resampledata(cfg, data);
end

%% Remove padding
cfg = [];
cfg.latency = [-params.pre params.post];
data = ft_selectdata(cfg, data); 

%% Convert to sensor definitions to cm
data.grad = ft_convert_units(data.grad,'cm');

%% Save bad trials
[~,idx]=ismember(data.sampleinfo,badtrl_jump,'rows');
badtrl_opm_jump = find(idx);
[~,idx]=ismember(data.sampleinfo,badtrl_std,'rows');
badtrl_opm_std = find(idx);
save(fullfile(save_path, [params.paradigm '_opm_badtrls']), ...
    'badtrl_opm_jump', ...
    'badtrl_opm_std',"-v7.3"); 

end

% for i = 1:138  
%     if sum(opm_raw.trial{1}(i,1:10)) == 0 & sum(opm_raw.trial{1}(i,1:200000))
%         disp(i)

%         % Create plot
%         fig = figure('Visible', 'off');  % Don’t show plot (faster)
%         plot(opm_raw.trial{1}(i,:));
%         title(sprintf('Sensor %d', i));

%         % Create filename with padded index (e.g., plot_001.png)
%         filename = sprintf('plot_%03d.png', i);

%         % Save the figure
%         saveas(fig, filename);

%         % Close the figure to save memory
%         close(fig);
%     end
% end
