function [ir1,ir2,ir3,x0,desired_point] = findnearestneighbour_scalar(irs,phi,delta,r,number_of_neighbours,X0)
    
% calculate the x-,y-,z-coordinates for the desired angles(phi,delta)
% and radius r
[x,y,z] = sph2cart(phi,delta,r);
desired_point = [x,y,z];


% calculate the x-,y-,z-coordinates for all known IRs
[x,y,z] = sph2cart(irs.apparent_azimuth(1,:), ...
irs.apparent_elevation(1,:),irs.distance(1,:));
x0 = [x;y;z];
    
if ~any(irs.apparent_azimuth(1,1)-irs.apparent_azimuth(1,2:end)) && ~any(irs.apparent_elevation(1,1)-irs.apparent_elevation(1,2:end))

    [~,idx] = findnearestneighbour_distance(x0,desired_point,X0);
    x0 = [x0(:,idx),[0;0;0]];
    ir1 = [irs.left(:,idx(1)) irs.right(:,idx(1))]; 
    ir2 = [irs.left(:,idx(2)) irs.right(:,idx(2))];
    ir3 = [0 0 0];
    
else

% calculate the scalarproduct between desired point and all points of 
% the HRIR set 
    scalarproduct = desired_point*x0;
    [~,idx] = sort(scalarproduct,'descend');
    idx = idx(1:number_of_neighbours);


% get the HRIR-positions corresponding to the maxima
    x0 = x0(:,idx);
    ir1 = [irs.left(:,idx(1)) irs.right(:,idx(1))];
    ir2 = [irs.left(:,idx(2)) irs.right(:,idx(2))];
    ir3 = [irs.left(:,idx(3)) irs.right(:,idx(3))]; 
    
end
fprintf('The indices of the found IRs with get_ir_new_version are: index1 = %d, index2 = %d\n',idx(1),idx(2));
fprintf('\n'); 

    
end
       



