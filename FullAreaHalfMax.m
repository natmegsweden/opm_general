function [peak_activation] = FullAreaHalfMax(sourcedistribution,sourcemodel,peak_latency,params)
%UNTITLED Calculates the full area at half max amplitude
%   Detailed explanation goes here

if length(peak_latency) == 2
    [~,i1] = min(abs(sourcedistribution.time-peak_latency(1)));
    [~,i2] = min(abs(sourcedistribution.time-peak_latency(2)));
    if iscell(sourcedistribution.avg.mom)
        dat = sourcedistribution.avg.pow(:,i1:i2);
    else
        dat = sourcedistribution.avg.mom(:,i1:i2).^2;
    end
    [~,i_latency] = max(max(dat,[],1)); % max of max across sources
    i_latency = i1-1+i_latency;
    latency = sourcedistribution.time(i_latency);
else
    [~,i_latency] = min(abs(sourcedistribution.time-peak_latency));
    latency = sourcedistribution.time(i_latency);
end

peak_activation = []; 
peak_activation.latency = latency;

if params.numdipoles == 2
    peak_activation.loc = nan(2,3);
    peak_activation.pow = nan(2,1);
    peak_activation.mom = nan(2,1);
    peak_activation.fahm = nan(2,1);
    peak_activation.halfmax_distribution = cell(2,1);
    peak_activation.target_region = nan(2,1);

    %% Left
    % Half max level
    brainstructures = find(startsWith(sourcemodel.brainstructurelabel,'L'));
    vertices = find(ismember(sourcemodel.brainstructure,brainstructures));
    if iscell(sourcedistribution.avg.mom)
        [peak_pow, i_maxsource] = max(abs(sourcedistribution.avg.pow(vertices,i_latency)));
        i_maxsource = vertices(i_maxsource);
        half_max = peak_pow/4;
        peak_loc = sourcemodel.pos(i_maxsource,:);
        peak_mom = norm(sourcedistribution.avg.mom{i_maxsource}(:,i_latency)); % power at max latency and source
    
        % Find triangles that have at least one point with amplitude >= half max
        halfmax_distribution = sourcedistribution.avg.pow(vertices,i_latency)>=half_max;
        i_halfmax_vertices = vertices(halfmax_distribution);
    else
        [peak_mom, i_maxsource] = max(abs(sourcedistribution.avg.mom(vertices,i_latency)));
        i_maxsource = vertices(i_maxsource);
        half_max = peak_mom/2;
        peak_loc = sourcemodel.pos(i_maxsource,:);
        peak_pow = peak_mom.^2; % power at max latency and source
    
        % Find triangles that have at least one point with amplitude >= half max   
        halfmax_distribution = abs(sourcedistribution.avg.mom(vertices,i_latency))>=half_max;
        i_halfmax_vertices = vertices(halfmax_distribution);
    end
    
    [triangles,~] = find(ismember(sourcemodel.tri,i_halfmax_vertices)); 
    triangles = sourcemodel.tri(triangles,:);
    
    % Sum area of triangles and divide by 3 (since its a triangle per point).
    fahm = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;  

    peak_activation.loc(1,:) = peak_loc;
    peak_activation.pow(1) = peak_pow;
    peak_activation.mom(1) = peak_mom;
    peak_activation.fahm(1) = fahm;
    peak_activation.halfmax_distribution{1} = i_halfmax_vertices;

    if isfield(params,'target_region')
        target_region = find(contains(sourcemodel.brainstructurelabel,params.target_region));
        peak_activation.target_region(1)  = sum(ismember(sourcemodel.brainstructure(i_halfmax_vertices),target_region))/length(i_halfmax_vertices);
    end

    %% Right
    brainstructures = find(startsWith(sourcemodel.brainstructurelabel,'R'));
    vertices = find(ismember(sourcemodel.brainstructure,brainstructures));
    if iscell(sourcedistribution.avg.mom)
        [peak_pow, i_maxsource] = max(abs(sourcedistribution.avg.pow(vertices,i_latency)));
        i_maxsource = vertices(i_maxsource);
        half_max = peak_pow/4;
        peak_loc = sourcemodel.pos(i_maxsource,:);
        peak_mom = norm(sourcedistribution.avg.mom{i_maxsource}(:,i_latency)); % power at max latency and source
    
        % Find triangles that have at least one point with amplitude >= half max
        halfmax_distribution = sourcedistribution.avg.pow(vertices,i_latency)>=half_max;
        i_halfmax_vertices = vertices(halfmax_distribution);
    else
        [peak_mom, i_maxsource] = max(abs(sourcedistribution.avg.mom(vertices,i_latency)));
        i_maxsource = vertices(i_maxsource);
        half_max = peak_mom/2;
        peak_loc = sourcemodel.pos(i_maxsource,:);
        peak_pow = peak_mom.^2; % power at max latency and source
    
        % Find triangles that have at least one point with amplitude >= half max   
        halfmax_distribution = abs(sourcedistribution.avg.mom(vertices,i_latency))>=half_max;
        i_halfmax_vertices = vertices(halfmax_distribution);
    end
    
    [triangles,~] = find(ismember(sourcemodel.tri,i_halfmax_vertices)); 
    triangles = sourcemodel.tri(triangles,:);
    
    % Sum area of triangles and divide by 3 (since its a triangle per point).
    fahm = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;  

    peak_activation.loc(2,:) = peak_loc;
    peak_activation.pow(2) = peak_pow;
    peak_activation.mom(2) = peak_mom;
    peak_activation.fahm(2) = fahm;
    peak_activation.halfmax_distribution{2} = i_halfmax_vertices;

    if isfield(params,'target_region')
        target_region = find(contains(sourcemodel.brainstructurelabel,params.target_region));
        peak_activation.target_region(2)  = sum(ismember(sourcemodel.brainstructure(i_halfmax_vertices),target_region))/length(i_halfmax_vertices);
    end

else
    % Half max level
    if iscell(sourcedistribution.avg.mom)
        [peak_pow, i_maxsource] = max(abs(sourcedistribution.avg.pow(:,i_latency)));
        half_max = peak_pow/4;
        peak_loc = sourcemodel.pos(i_maxsource,:);
        peak_mom = norm(sourcedistribution.avg.mom{i_maxsource}(:,i_latency)); % power at max latency and source
    
        % Find triangles that have at least one point with amplitude >= half max
        halfmax_distribution = sourcedistribution.avg.pow(:,i_latency)>=half_max;
        i_halfmax_vertices = find(halfmax_distribution);
    else
        [peak_mom, i_maxsource] = max(abs(sourcedistribution.avg.mom(:,i_latency)));
        half_max = peak_mom/2;
        peak_loc = sourcemodel.pos(i_maxsource,:);
        peak_pow = peak_mom.^2; % power at max latency and source
    
        % Find triangles that have at least one point with amplitude >= half max   
        halfmax_distribution = abs(sourcedistribution.avg.mom(:,i_latency))>=half_max;
        i_halfmax_vertices = find(halfmax_distribution);
    end
    
    [triangles,~] = find(ismember(sourcemodel.tri,i_halfmax_vertices)); 
    triangles = sourcemodel.tri(triangles,:);
    
    % Sum area of triangles and divide by 3 (since its a triangle per point).
    fahm = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;  
    
    peak_activation.loc = peak_loc;
    peak_activation.pow = peak_pow;
    peak_activation.mom = peak_mom;
    peak_activation.fahm = fahm;
    peak_activation.halfmax_distribution = i_halfmax_vertices;
    
    if isfield(params,'target_region')
        target_region = find(contains(sourcemodel.brainstructurelabel,params.target_region));
        peak_activation.target_region  = sum(ismember(sourcemodel.brainstructure(i_halfmax_vertices),target_region))/length(i_halfmax_vertices);
    end
end
end