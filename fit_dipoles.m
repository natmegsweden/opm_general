function [dipole] = fit_dipoles(save_path,timelocked,headmodel,mri, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numdipoles = params.numdipoles;

colors = [[0 0.4470 0.7410]; % blue
    [0.8500 0.3250 0.0980]; % red
    [0.9290 0.6940 0.1250]; % yellow
    [0.4940 0.1840 0.5560]; % purple
    [0.4660 0.6740 0.1880]; % green
    [0.6350 0.0780 0.1840]]; % light blue

cfg              = [];
cfg.resolution   = 0.5;
cfg.tight        = 'yes';
cfg.inwardshift  = 0;
cfg.headmodel    = headmodel;
sourcemodel    = ft_prepare_sourcemodel(cfg);

%% Fit dipoles
% MEG
cfg = [];
cfg.gridsearch      = 'yes';           
cfg.numdipoles      = 1;                
cfg.sourcemodel     = sourcemodel;           
cfg.headmodel       = headmodel;    
cfg.senstype        = 'meg';            
cfg.channel         = params.chs;  
cfg.numdipoles      = numdipoles;
if numdipoles == 2
    cfg.symmetry    = 'x';
end
cfg.nonlinear       = 'yes';           
cfg.latency         = timelocked.peak_latency;% + [-0.005 0.005];
cfg.dipfit.checkinside = 'yes';
dipole = ft_dipolefitting(cfg, timelocked);
close all

%% Make sure symmetric dipoles are ordered the same
if numdipoles == 2
    if dipole.dip.pos(1,1) > 0 % first dipole is on the right -> flip order
        dipole.dip.pos([2 1],:) = dipole.dip.pos([1 2],:);
        dipole.dip.mom([4 5 6 1 2 3]) = dipole.dip.mom([1 2 3 4 5 6]);
    end
end

%% Plot triggers jointly
pos = zeros(numdipoles,3);
ori = zeros(numdipoles,3);
mom = zeros(numdipoles,1);
rv = zeros(numdipoles,1);
for i_dip = 1:numdipoles
    pos(i_dip,:) =dipole.dip.pos(i_dip,:);
    [mom(i_dip,1),idx] = max(vecnorm(dipole.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
    ori(i_dip,:) = dipole.dip.mom((1:3)+((i_dip-1)*3),idx);
    rv(i_dip,:) = dipole.dip.rv;
end
dipole.dip.peak_mom = mom;
mean_pos = mean(pos,1);

tmp = pos;
tmp(:,1) = mean_pos(1)-1;

h=figure;
h.Position = [100 100 800 800];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i_trigger = 1:numdipoles
    ft_plot_dipole(tmp(i_trigger,:),ori(i_trigger,:),'color',colors(i_trigger,:))
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i_trigger = 1:numdipoles
    ft_plot_dipole(tmp(i_trigger,:),ori(i_trigger,:),'color',colors(i_trigger,:))
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i_trigger = 1:numdipoles
    ft_plot_dipole(tmp(i_trigger,:),ori(i_trigger,:),'color',colors(i_trigger,:))
end
hold off
axis equal; axis tight; axis vis3d; axis off
h.Position = [100 100 800 900];
subplot(2,2,4);
axis off; % Turn off axis
text(0, 0.88, [params.paradigm ' - ' params.modality], 'FontWeight', 'bold');
text(0, 0.8, 'latency (ms): ');
text(0.3, 0.8, num2str(timelocked.peak_latency*1e3));
text(0, 0.75, 'dipole: ');
text(0, 0.7, 'pos (mm): '); 
text(0, 0.65, 'mom (nAm): '); 
text(0, 0.6, 'rv: '); 
for i_dip = 1:numdipoles
    text(0.3 + (i_dip-1)*0.35, 0.75, num2str(i_dip)); 
    text(0.3 + (i_dip-1)*0.35, 0.7, [num2str(10*pos(i_dip,1),'%.1f') ' / ' num2str(10*pos(i_dip,2),'%.1f') ' / ' num2str(10*pos(i_dip,3),'%.1f')]); 
    text(0.3 + (i_dip-1)*0.35, 0.65, num2str(1e9*1e-4*mom(i_dip),'%.3f')); 
    text(0.3 + (i_dip-1)*0.35, 0.6, num2str(rv(i_dip),'%.3f')); 
end
saveas(h, fullfile(save_path, 'figs', [params.paradigm '_' params.modality '_dipfit_mri.jpg']))
close all

%% Save 
save(fullfile(save_path, [params.paradigm '_' params.modality '_dipoles']), 'dipole'); disp('done');

end