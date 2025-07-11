function fit_mne(save_path, timelocked, headmodel, sourcemodel, sourcemodel_inflated, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(params,'plot_inflated')
    params.plot_inflated = false;
end

%% Prepare leadfields
headmodel = ft_convert_units(headmodel,'cm');
sourcemodel = ft_convert_units(sourcemodel,'cm');
if params.source_fixedori
    sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
end

%% MNE invserse
mne_peak = cell(length(params.trigger_code),length(params.peaks));

%% Leadfield
cfg = [];
cfg.grad             = timelocked{1}.grad; % sensor positions
cfg.channel          = params.chs;
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
%cfg.normalize        = 'yes';
leadfield = ft_prepare_leadfield(cfg,timelocked{1});

for i_trigger = 1:length(params.trigger_code)
    params.i_trigger = i_trigger;
    cov = '';
    if isfield(params,'use_cov') && strcmp(params.use_cov,'all')
        timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_all;
        cov = '_covAll';
    elseif isfield(params,'use_cov') && strcmp(params.use_cov,'resting_state')
        cov = '_covRS';
        if size(timelocked{i_trigger}.cov_RS,1) >= size(timelocked{i_trigger}.cov,1)
            timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_RS;
        else
            disp([params.sub '_' params.trigger_labels{params.i_trigger} ': ' num2str(size(timelocked{i_trigger}.cov_RS,1)) ' vs ' num2str(size(timelocked{i_trigger}.cov,1))])
            peak{i_trigger,1} = [params.sub '_' params.trigger_labels{params.i_trigger} ': ' num2str(size(timelocked{i_trigger}.cov_RS,1)) ' vs ' num2str(size(timelocked{i_trigger}.cov,1))];
            continue
        end
    elseif isfield(params,'use_cov') && strcmp(params.use_cov,'empty_room')
        cov = '_covER';
        if ~isfield(timelocked{i_trigger},'cov_ER')
            continue % skip if no empty room covariance available
        end
        timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_ER;
        if size(timelocked{i_trigger}.cov_ER) == size(timelocked{i_trigger}.cov)
            timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_ER;
        else
            peak{i_trigger,1} = [params.sub '_' params.trigger_labels{params.i_trigger} ': ' num2str(size(timelocked{i_trigger}.cov_ER,1)) ' vs ' num2str(size(timelocked{i_trigger}.cov,1))];
            continue
        end
    end

    %% Source reconstruction
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = params.chs;
    tmp = ft_sourceanalysis(cfg, timelocked{i_trigger});
    tmp.tri = sourcemodel.tri;

    if iscell(tmp.avg.mom) && size(tmp.avg.mom{1},1)~=3
        if size(tmp.avg.mom,1) == 1
            tmp.avg.mom = cell2mat(tmp.avg.mom');
        else
            tmp.avg.mom = cell2mat(tmp.avg.mom);
        end
    end
    
    h = figure;
    plot(tmp.time*1e3,tmp.avg.pow)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.trigger_labels{params.i_trigger}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_mne_sourcepow_ph-' params.trigger_labels{params.i_trigger} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
       peak{i_trigger,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = peak{i_trigger,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['SQUID-MAG (FAHM=' num2str(peak{i_trigger,i_peak}.fahm,3) 'cm^2; t=' num2str(round(peak{i_trigger,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak{i_trigger,i_peak}.label '_mne_ph' params.trigger_labels{i_trigger} cov '.jpg']))
        close all
    end

    clear tmp

end
save(fullfile(save_path, [params.modality '_mne_peaks' cov]), 'peak'); 

end