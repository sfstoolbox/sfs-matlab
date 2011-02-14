function kemarmit2mat(irspath)
%KEMARMIT2MAT converts IRs data given by KEMAR MIT to our mat-file based format
%
%   Usage: kemarmit2mat(irspath);
%
%   Input options:
%       irspath - path to the directory containing the IR data
%
%   KEMARMIT2MAT(irspath) converts the IRs data given by the KEMAR MIT
%   measurement and stored at the given irspath in our own mat-file based 
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
    ['MIT anechoic measurements with KEMAR',...
     'Used elevation angle: 0°; azimuth resolution: 5°. ', ...
     'For further information have a look at: ',...
     'http://sound.media.mit.edu/resources/KEMAR.html'];
irs.head = 'KEMAR DB-4004';
irs.room = 'Anechoic chamber of MIT';
irs.loudspeaker = 'Optimus Pro 7';
irs.distance = 1.7;
irs.fs = 44100;
irs.head_position = [0 0 0];
irs.head_reference = [0 1 0];
irs.source_position = [0 1.7 0];
irs.source_reference = [0 0 0];
irs.head_azimuth = NaN;
irs.head_elevation = NaN;

% Read the data
elevs = -40:10:80;
steps = [6.43,6,5,5,5,5,5,6,6.43,8,10,15,30];
idx = 1;
for jj = 1:length(elevs)
    for phi = 0:steps(jj):180
        irfile = sprintf('%s/elev%i/H%ie%03.0fa.wav',irspath,elevs(jj),...
            elevs(jj),phi);
        ir = wavread(irfile);
        irs.apparent_azimuth(idx) = correct_azimuth(rad(-phi));
        irs.apparent_elevation(idx) = correct_elevation(rad(elevs(jj)));
        irs.head_azimuth(idx) = correct_azimuth(rad(phi));
        irs.head_elevation(idx) = correct_elevation(rad(-elevs(jj)));
        irs.left(:,idx) = ir(:,1);
        irs.right(:,idx) = ir(:,2);
        idx = idx+1;
    end
    for phi = steps(jj):steps(jj):180-steps(jj)
        irfile = sprintf('%s/elev%i/H%ie%03.0fa.wav',irspath,elevs(jj),...
            elevs(jj),phi);
        ir = wavread(irfile);
        irs.apparent_azimuth(idx) = correct_azimuth(rad(phi));
        irs.apparent_elevation(idx) = correct_elevation(rad(elevs(jj)));
        irs.head_azimuth(idx) = correct_azimuth(rad(-phi));
        irs.head_elevation(idx) = correct_elevation(rad(-elevs(jj)));
        irs.left(:,idx) = ir(:,2);
        irs.right(:,idx) = ir(:,1);
        idx = idx+1;
    end
end
irfile = sprintf('%s/elev90/H90e000a.wav',irspath);
ir = wavread(irfile);
irs.apparent_azimuth(idx) = correct_azimuth(0);
irs.apparent_elevation(idx) = correct_elevation(rad(90));
irs.head_azimuth(idx) = correct_azimuth(0);
irs.head_elevation(idx) = correct_elevation(rad(-90));
irs.left(:,idx) = ir(:,1);
irs.right(:,idx) = ir(:,2);

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
outfile = sprintf('%s/KEMAR_MIT.mat',outdir);
save('-v7',outfile,'irs');
