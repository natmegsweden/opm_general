function h = plot_source_distribution(srcdist, peak, params)
%plot_source_distribution Function for plotting source distributions
%   3D plot of a source distribution on a sourcemodel. Expects peak to
%   contain the fields peak.latency (latency at which to display source
%   activation) and peak.fahm.
%   Use params.mne_view to select how the figure should be displayed. 
%   params.mne_view options:
%   [az el]     Azimutal and elevation angles determining the view angle.
%   'sides'     Display sourcemodel from left and right to show both
%               hemisphres
%   default     If peak structure from FullAreaHalfMax function exists,
%               orient camera to focus on peak activation
if isfield(params,'mne_view') && isnumeric(params.mne_view) && length(params.mne_view) == 2
    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = peak.latency;
    h = figure;
    h.Position(3) = round(h.Position(3)*1.2);
    ft_sourceplot(cfg, srcdist)
    lighting gouraud
    material dull
    view(params.mne_view(1),params.mne_view(2))
    if length(peak.fahm)>1
        title([params.paradigm '-' params.modality ' (FAHM_R=' num2str(peak.fahm(2),3) 'cm^2; FAHM_L=' num2str(peak.fahm(1),3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    else
        title([params.paradigm '-' params.modality ' (FAHM=' num2str(peak.fahm,3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    end
elseif isfield(params,'mne_view') && strcmp(params.mne_view,'sides')
    viewangles = [90 0 -90 0];
    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = peak.latency;
    h = figure;
    h.Position(3) = round(h.Position(3)*1.6);
    subplot(1,2,1); % right hemisphere
    cfg.figure = h;
    ft_sourceplot(cfg, srcdist)
    material dull
    view(viewangles(1),viewangles(2))
    camlight();
    lighting gouraud
    if length(peak.fahm)>1
        title([params.paradigm '-' params.modality ' (FAHM_R=' num2str(peak.fahm(2),3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    else
        title([params.paradigm '-' params.modality ' (FAHM=' num2str(peak.fahm,3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    end
    subplot(1,2,2); % left hemisphere
    cfg.figure = h;
    ft_sourceplot(cfg, srcdist)
    material dull
    view(viewangles(3),viewangles(4))
    camlight()
    lighting gouraud
    if length(peak.fahm)>1
        title([params.paradigm '-' params.modality ' (FAHM_L=' num2str(peak.fahm(1),3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    else
        title([params.paradigm '-' params.modality ' (FAHM=' num2str(peak.fahm,3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    end
else
    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = peak.latency;
    h = figure;
    h.Position(3) = round(h.Position(3)*1.2);
    ft_sourceplot(cfg, srcdist)
    lighting gouraud
    material dull
    if isfield(peak.loc)
        center = mean(srcdist.pos,1);
        direction = mean(peak.loc,1) - center;
        direction = direction / norm(direction);
        distance = 200;
        campos(center + distance * direction);
        camtarget(peak.loc);
        camup([0 0 1]);
        camproj('perspective');
    end
    if length(peak.fahm)>1
        title([params.paradigm '-' params.modality ' (FAHM_R=' num2str(peak.fahm(2),3) 'cm^2; FAHM_L=' num2str(peak.fahm(1),3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    else
        title([params.paradigm '-' params.modality ' (FAHM=' num2str(peak.fahm,3) 'cm^2; t=' num2str(round(peak.latency*1e3)) 'ms)'])
    end
end
end