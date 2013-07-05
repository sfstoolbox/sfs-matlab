function fabian3D2mat(irsset)
%FABIAN3D2mat converts 3D IRs data given by the FABIAN HRIR measurement  
% from Oldenburg to our mat-file based format
%
%   Usage: fabian3D2mat(irsset) with irsset = load('az0.mat');
%
%   Note: irsset has to be an *.mat file containing an irs-struct 
% 
%   Input options:
%       irsset  - IR sets measured with FABIAN. Currently the 
%                 following are available:
%					'az0.mat'  
%
%   FABIAN3D2mat(irsset) converts the IRs data given by the irsset
%   in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format
%
% AUTHOR: Vincent Kuscher


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
outdir = 'ir_databases';


%% ===== Computation =====================================================
irs = irsset;
irs.apparent_azimuth = correct_azimuth(rad(irs.source_azimuth)');
irs.apparent_elevation = correct_elevation(rad(irs.source_elevation)');
irs.distance = irs.distance';
irs = rmfield(irs,'source_azimuth');
irs = rmfield(irs,'source_elevation');
irs.ears = 'ear';
irs.head_position = [0;0;0];
irs.head_reference = [0;1;0];

irs.source_reference = [0;0;0];
irs.head_elevation = NaN;
irs.torso_azimuth = NaN;
irs.torso_elevation = NaN;
[x,y,z] =sph2cart(irs.apparent_azimuth,irs.apparent_elevation,irs.distance);
irs.source_position = [x;y;z];

irs = rmfield(irs,'microphone');
irs = rmfield(irs,'itd');

irs = order_irs_fields(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
OutputName = 'FABIAN_3D_anechoic';
save([OutputName '.mat'], 'irs');
