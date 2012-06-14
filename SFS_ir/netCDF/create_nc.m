% Create a netCDF HRTF file
% see: http://www.unidata.ucar.edu/software/netcdf/

% Read irs
addirspath;
irs = read_irs('QU_KEMAR_anechoic_3m.mat');

% dimensions
samples = length(irs.left);
points = length(irs.apparent_azimuth);
channels = 2;
positions = 6;


%% === Data ===
% ir
data(:,:,1) = irs.left;
data(:,:,2) = irs.right;
% % receiver
% receiver = zeros(position,points,channels);
% receiver(:,:,1) = repmat([-0.08 0 0 0 0]',[1 points]);
% % head
% head = zeros(position,points);
% if length(irs.head_azimuth)>1 && length(irs.head_elevation)>1
%     for ii = 1:points
%         head(:,ii) = [0 0 0 irs.head_azimuth(ii) irs.head_elevation(ii)];
%     end
% elseif length(irs.head_azimuth)>1
%     for ii = 1:points
%         head(:,ii) = [0 0 0 irs.head_azimuth(ii) 0];
%     end
% elseif length(irs.head_elevation)>1
%     for ii = 1:points
%         head(:,ii) = [0 0 0 0 irs.head_elevation(ii)];
%     end
% else
%     for ii = 1:points
%         head(:,ii) = [0 0 0 0 0];
%     end
% end
% body
receiver = zeros(positions,points,channels);
if length(irs.torso_azimuth)>1 && length(irs.torso_elevation)>1
    for ii = 1:points
        [nx ny nz] = ...
            sph2cart(irs.torso_azimuth(ii)+pi/2,irs.torso_elevation(ii),1);
        receiver(:,ii,:) = repmat([0 0 0 nx ny nz],[2 1])';
    end
elseif length(irs.torso_azimuth)>1
    for ii = 1:points
        [nx ny nz] = sph2cart(irs.torso_azimuth(ii)+pi/2,0,1);
        receiver(:,ii,:) = repmat([0 0 0 nx ny nz],[2 1])';
    end
elseif length(irs.torso_elevation)>1
    for ii = 1:points
        [nx ny nz] = sph2cart(0,irs.torso_elevation(ii),1);
        receiver(:,ii) = repmat([0 0 0 nx ny nz],[2 1])';
    end
else
    receiver(:,:,1) = repmat([0 0 0 0 1 0]',[1 points]);
    receiver(:,:,2) = repmat([0 0 0 0 1 0]',[1 points]);
end
% source
source = zeros(positions,points);
source(:,:) = repmat([irs.source_position' 0 -1 0]',[1 points]);


%% === Save file ===
fid = 'QU_KEMAR_anechoic_3m.nc';

if exist(fid,'file')
    delete(fid);
end

% data
nccreate(fid, ...
         '/data', ...
         'Dimensions',{'samples',samples, ...
                       'points',points, ...
                       'channels',channels}, ...
         'DataType','double', ...
         'Format','netcdf4', ...
         'DeflateLevel',9);
ncwrite(fid,'/data',data);
ncwriteatt(fid,'/data','Type','FIR');
ncwriteatt(fid,'/data','SamplingRate',44100);
% receiver
nccreate(fid, ...
         '/receiver', ...
         'Dimensions',{'positions',positions, ...
                       'points',points, ...
                       'channels',channels}, ...
         'DataType','double', ...
         'Format','netcdf4', ...
         'DeflateLevel',9);
ncwrite(fid,'/receiver',receiver);
ncwriteatt(fid,'/receiver','Description',[irs.head,', ears: ',irs.ears]);
% % head
% nccreate(fid, ...
%          '/head', ...
%          'Dimensions',{'position',position, ...
%                        'points',points}, ...
%          'DataType','double', ...
%          'Format','netcdf4', ...
%          'DeflateLevel',9);
% ncwrite(fid,'/head',head);
% ncwriteatt(fid,'/head','Description',irs.head);
% % body
% nccreate(fid, ...
%          '/body', ...
%          'Dimensions',{'position',position, ...
%                        'points',points}, ...
%          'DataType','double', ...
%          'Format','netcdf4', ...
%          'DeflateLevel',9);
% ncwrite(fid,'/body',body);
% ncwriteatt(fid,'/body','Description',irs.head);
% source
nccreate(fid, ...
         '/source', ...
         'Dimensions',{'positions',positions, ...
                       'points',points}, ...
         'DataType','double', ...
         'Format','netcdf4', ...
         'DeflateLevel',9);
ncwrite(fid,'/source',source);
ncwriteatt(fid,'/source','Description',irs.source);
% metadata
ncwriteatt(fid,'/','OrientationType','navigational');
ncwriteatt(fid,'/','PositionType','cartesian');
ncwriteatt(fid,'/','Description',irs.description);
ncwriteatt(fid,'/','SubjectID',irs.head);
ncwriteatt(fid,'/','AuthorName','H. Wierstorf, M. Geier, S. Spors');
ncwriteatt(fid,'/','AuthorContact','hagen.wierstorf@tu-berlin.de');
ncwriteatt(fid,'/','License','Creative Commons Attribution-ShareAlike 3.0');
ncwriteatt(fid,'/','Organization','Quality and Usability Lab, T-Labs, TU Berlin');
ncwriteatt(fid,'/','DatabaseName','QU_KEMAR_anechoic_3m');
ncwriteatt(fid,'/','RoomDescription',irs.room);

 