function wittek2mat(irspath)
%WITTEK2MAT converts IRs data given by Wittek to our mat-file based format
%
%   Usage: wittek2mat(irspath);
%
%   Input options:
%       irspath - path to the directory containing the IR data
%
%   WITTEK2MAT(irspath) converts the IRs data given by the measurements from
%   Wittek and stored at the given irspath in our own mat-file based 
%   format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargdir(irspath);


%% ===== Computation =====================================================
outdir = 'ir_databases';
% Initialize a new IR struct
irs = new_irs();
irs.description = ...
    ['Robert Wittek measurements with KEMAR in a Studio.',...
     'Used elevation angle: 0 deg; azimuth resolution: 1 deg. '];
irs.head = 'KEMAR';
irs.room = 'Studio';
irs.loudspeaker = 'dummy';
irs.distance = 1;
irs.fs = 44100;
irs.head_position = [0 0 0];
irs.head_reference = [0 1 0];
irs.source_position = [0 1 0];
irs.source_reference = [0 0 0];
irs.head_azimuth = NaN;
irs.head_elevation = NaN;
irs.torso_elevation = NaN;

% Read the data
warning('off', 'all');
idx = 1;
for phi = 0:1:359

    irfile = sprintf('%s/Studio360_KHneu_L3600_%i.mat',irspath,phi*10);
    load(irfile);
    irs.apparent_azimuth(idx) = correct_azimuth(rad(phi));
    irs.apparent_elevation(idx) = correct_elevation(0);
    irs.torso_azimuth(idx) = -correct_azimuth(rad(phi));
    % FIXME: hrtf is of class data. How do I handle this?
    irs.left(:,idx) = hrtf.data(:,1);
    irs.right(:,idx) = hrtf.data(:,2);
    idx = idx+1;

end

% Reorder fields
irs = order_irs_fields(irs);
% Reorder entries
irs = correct_irs_angle_order(irs);
% Check irs format
check_irs(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
outfile = sprintf('%s/KEMAR_Wittek.mat',outdir);
save('-v7',outfile,'irs');
