function [inside] = insidePointcloud(points, point_cloud)
%UNTITLED3 Checks whether point is inside a closed convex mesh
%   point = 1x3
%   pos : Nx3
n_pos = size(point_cloud,1);
n_points = size(points,1);
inside = false(n_points,1);
% find triangle centers
centroid = mean(point_cloud); 

% find closest point from pos
for i = 1:size(points,1)
    [~, i_min] = min(vecnorm(point_cloud-repmat(points(i,:),[n_pos 1]),2,2));
    if dot(centroid-point_cloud(i_min,:),points(i,:)-point_cloud(i_min,:))>0
        inside(i) = true;
    end
end
end