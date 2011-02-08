function bkoldenburg2mat(irsset,irspath)
%BKOLDENBURG2MAT converts IRs from BK dummy head from Oldenburg to irs format
%   Usage: bkoldenburg2mat(irspath);
%
%   Input options:
%       irsset - IR sets:
%                   'RAR_08m'   - HRIR of RAR with 0.8 m distance.
%                   'RAR_3m'    - HRIR of RAR with 3 m distance.
%       irspath - path to the directory containing the IR data
%
%   BKOLDENBURG2MAT(irsset,irspath) converts the IRs data given by the irsset
%   and measured with a B&K dummy head in Oldenburg and stored at the given 
%   irspath in our own mat-file based format. For format details have a look 
%   at the IR_format.txt file.

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargchar(irsset);
isargdir(irspath);


%% ===== Computation =====================================================
outdir = 'ir_databases';

% Initialize a new IR struct
irs = new_irs();
% Add common entries
irs.fs = 48000;
irs.loudspeaker = 'Tannoy 800A LH';
irs.room = 'Anechoic chamber of Carl-von-Ossietzky Universität Oldenburg';
irs.head = 'Brüel & Kjaer Type 4128C';
irs.head_elevation = NaN;
irs.head_azimuth = NaN;
irs.head_position = [0 0 0];
irs.head_reference = [0 1 0];
irs.source_reference = [0 0 0];

if strcmp(irsset,'RAR_08m')
    irs.description = ...
        ['Anechoic measurements with a B&K dummy head in the anechoic ',...
         'chamber of Carl-von-Ossietzky Universität Oldenburg. ',...
         'Used elevation angle: -10 degree,0 degree, 10 degree ,20 degree; ',...
         'azimuth resolution: 5 degree. Distance: 0.8 m. ', ...
         'For further information have a look at: ',...
         'http://medi.uni-oldenburg.de/hrir/index.html'];
    irs.distance = 0.8;
    irs.source_position = [0 0.8 0];
elseif strcmp(irsset,'RAR_3m')
    irs.description = ...
        ['Anechoic measurements with a B&K dummy head in the anechoic ',...
         'chamber of Carl-von-Ossietzky Universität Oldenburg. ',...
         'Used elevation angle: -10 degree,0 degree, 10 degree ,20 degree; ',...
         'azimuth resolution: 5 degree. Distance: 3 m. ', ...
         'For further information have a look at: ',...
         'http://medi.uni-oldenburg.de/hrir/index.html'];
    irs.distance = 3;
    irs.source_position = [0 3 0];
end

% Read data
idx = 1;
for delta = -10:10:20
    for azimuth = -180:5:175
        irfile = sprintf('%s/anechoic_distcm_%i_el_%i_az_%i.wav',...
            irspath,irs.distance*100,delta,azimuth);
        sig = wavread(irfile);
        irs.left(:,idx) = sig(:,1);
        irs.right(:,idx) = sig(:,2);
        irs.apparent_azimuth(idx) = -correct_azimuth(rad(azimuth));
        irs.torso_azimuth(idx) = correct_azimuth(rad(azimuth));
        irs.apparent_elevation(idx) = correct_elevation(rad(delta));
        irs.torso_elevation(idx) = -correct_elevation(rad(delta));
        idx = idx+1;
    end
end


% Reorder entries
irs = correct_irs_angle_order(irs);
irs = order_irs_fields(irs);
check_irs(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
outfile = sprintf('%s/BKOldenburg_%s.mat',outdir,irsset);
save('-v7',outfile,'irs');
