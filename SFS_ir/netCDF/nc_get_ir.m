function [ir,metadata] = nc_get_ir(file,angle)
%NC_GET_IR returns the IR for the given angle from file
%   Usage: ir = nc_get_ir(file,angle)
%

% AUTHOR: Hagen Wierstorf

ir = [];
metadata = struct;

if isoctave %% === Octave ===

    % NOTE: in octave the order of dimension within the variables is converted!

    % read position vectors
    nc = netcdf(file,'r');
    receiver = nc{'receiver'};
    metadata.receiver = nc{'receiver'}.Description;
    source = nc{'source'};
    metadata.source = nc{'source'}.Description;

    % position of head
    head_pos = head(:,1:3) + body(:,1:3);
    src_pos = source(:,1:3);
    % calculate apparent azimuth
    tmp = head_pos-src_pos;
    [src_phi,src_th,distance] = cart2sph(tmp(:,1),tmp(:,2),tmp(:,3));
    head_phi = head(:,4)+body(:,4);
    apparent_azimuth = src_phi + head_phi + pi/2;
    degree(apparent_azimuth)

    idx = find(apparent_azimuth>=angle,1,'first');

    % read corresponding IR
    metadata.distance = distance(idx);
    ir = nc{'ir'}(:,idx,:);
    ir = squeeze(ir);
    ir = ir';


else %% === Matlab ===

    % read position vectors
    receiver = ncread(file,'/receiver');
    metadata.receiver = ncreadatt(file,'/receiver','Description');
    %head = ncread(file,'/head');
    %metadata.head = ncreadatt(file,'/head','Description');
    %body = ncread(file,'/body');
    %metadata.body = ncreadatt(file,'/body','Description');
    source = ncread(file,'/source');
    metadata.source = ncreadatt(file,'/source','Description');

    % position of head
    head_pos = receiver(1:3,:,1);
    src_pos = source(1:3,:);
    % calculate apparent azimuth
    tmp = head_pos-src_pos;
    [src_phi,~,distance] = cart2sph(tmp(1,:),tmp(2,:),tmp(3,:));
    nx = receiver(4,:,1);
    ny = receiver(5,:,1);
    nz = receiver(6,:,1);
    head_phi = cart2sph(nx,ny,nz);
    apparent_azimuth = src_phi + head_phi;

    idx = find(apparent_azimuth>=angle,1,'first');

    % Get corresponding IR
    metadata.distance = distance(idx);
    ir = ncread(file,'/data',[1 idx 1],[Inf 1 Inf]);
    ir = squeeze(ir);

end
