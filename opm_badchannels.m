function [badchs, badchs_flat, badchs_std, badchs_neighbors, badchs_outlier] = opm_badchannels(data, trl, params)
%opm_badchannels Detects channels that are flat, have low correlation with
%thei_chs_gradeighbors or show a lot of jumping artifacts.
%   cfg.z_threshold
%   cfg.corr_threshold
%   cfg.n_neighbors
%   cfg.njump_threshold

std_threshold = ft_getopt(params, 'std_threshold', 5e-9);
n_neighbors     = ft_getopt(params, 'n_neighbors', 4);
corr_threshold  = ft_getopt(params, 'corr_threshold', 0.6);
z_threshold     = ft_getopt(params, 'z_threshold', 20);
njump_threshold = ft_getopt(params, 'njump_threshold', 0.05);

cfg = [];
cfg.channel = '*bz*';
data = ft_selectdata(cfg, data);

chs = find(contains(data.label,'_bz*'));

%% Find channels with flat segments or high std
cfg = [];
cfg.length = 1;
data_seg = ft_redefinetrial(cfg,data);
for i_trl = 1:length(data_seg.trial)
    trl_std(:,i_trl) = std(data_seg.trial{i_trl},0,2);
end
badchs_flat = find(any(trl_std<1e-15,2));
badchs_std = find(mean(trl_std,2)>std_threshold);

%% Neighbors
goodchs = setdiff(chs,[badchs_flat; badchs_std]);

cfg = [];
cfg.resamplefs = 200;
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.lpinstabilityfix  = 'reduce';
data_lp = ft_resampledata(cfg,data);

cfg = [];
cfg.length = 1;
data_lp = ft_redefinetrial(cfg,data_lp);

% Create neighbor structure
n_chs = length(data_lp.label);
chanpos = data_lp.grad.chanpos;
neighbors = zeros(n_chs,n_neighbors);
for i = 1:size(chanpos,1)
        [~,tmp]= sort(vecnorm(chanpos(goodchs,:)-repmat(chanpos(i,:),[length(goodchs) 1]),2,2));
        neighbors(i,:) = goodchs(tmp(2:(n_neighbors+1)));
end
neighborscorr = zeros(n_chs,n_neighbors,length(data_lp.trial));
for i_trl = 1:length(data_lp.trial)
    dat = data_lp.trial{i_trl};
    for i = 1:n_chs
        for j = 1:n_neighbors
                tmp2 = corrcoef(dat(i,:),dat(int32(neighbors(i,j)),:));
                neighborscorr(i,j,i_trl) = abs(tmp2(1,2));
        end
    end 
end
badchs_neighbors = find(max(mean(neighborscorr,3),[],2)<corr_threshold); % bad if no neighbors exceed correlation threshold
badchs = [badchs_flat; badchs_std; badchs_neighbors];

%% Epoch
cfg = [];
cfg.trl = trl;
data_epo = ft_redefinetrial(cfg,data);
cfg = [];
cfg.demean = 'yes';
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = [40 50 60 80 100 ];
data_epo = ft_preprocessing(cfg,data_epo);

%% Spectrum
goodchs = setdiff(chs,badchs);

cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [3 70];
freq = ft_freqanalysis(cfg, data_epo);
freq.powspctrm = sqrt(freq.powspctrm);
threshold = mean(freq.powspctrm(goodchs,:),1)+3*std(freq.powspctrm(goodchs,:),0,1);
badchs_outlier = find(mean(freq.powspctrm>repmat(threshold,[size(freq.powspctrm,1) 1]),2)>0.5); % more than half the frequencies above threshold

badchs = [badchs_flat; badchs_std; badchs_neighbors; badchs_outlier];

% Convert to channel labels
badchs = data.label(chs(badchs));
badchs_flat = data.label(chs(badchs_flat));
badchs_std = data.label(chs(badchs_std));
badchs_neighbors = data.label(chs(badchs_neighbors));
badchs_outlier = data.label(chs(badchs_outlier));