function aachen2mat(irsset)

%AACHEN2MAT converts IRs data given by the Aachen_Impulse_Response_Database_v1.2 to our mat-file based format
%
%   Usage: aachen2mat(irsset);
%
%   Input options:
%       irsset  - IR sets measured with VariSphear and KEMAR. Currently the 
%                 following are available:
%					'stairway_1m' - Binaural RIR of Stairways with 1m distance
%					'stairway_2m' - Binaural RIR of Stairway with 1m distance
%					'stairway_3m' - Binaural RIR of Stairways with 1m distance
%
%                 NOTE: you still have to give the matching path to the given
%                 data set!
%
%       irspath - path to the directory containing the IR data
%
%   AACHEN2MAT(irsset,irspath) converts the IRs data given by the irsset
%   in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format
%
% AUTHOR: Lars-Erik Riechert

%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargchar(irsset);

%% ===== Computation =====================================================
% Dir to save the IRs
outdir = 'ir_databases';

% Initialize a new IR struct
irs = new_irs();

% Add common struct entries
irs.fs = 48000;
irs.source = 'NaN';
irs.head_elevation = NaN;
irs.torso_elevation = NaN;
irs.head_position = [0 0 0]';
irs.head_reference = [0 1 0]';
irs.source_reference = [0 0 0]';

irs.description = ...
        ['Aachen_Impulse_Response_Database_v1.2; azimuth resolution: 15°. '];
irs.ears = 'NaN';
irs.torso_azimuth = NaN;
irs.room_corners = [0 0 -1.7;]';
	
% Setup for the load_air script 
airpar.fs = 48e3;
airpar.rir_type = 1;
airpar.room = 5;




if strcmp(irsset,'stairway_1m')
	irs.source_position = [0 1 0]';
	airpar.rir_no = 1;
	irs.distance = 1;
	airpar.head = 1;
	irs.head = 'some head';
elseif strcmp(irsset,'stairway_2m')
	irs.source_position = [0 2 0]';
	airpar.rir_no = 2;
	irs.distance = 2;
	airpar.head = 1;
	irs.head = 'some head';
elseif strcmp(irsset,'stairway_3m')
	irs.source_position = [0 3 0]';
	airpar.rir_no = 3;
	irs.distance = 3;
	airpar.head = 1;
	irs.head = 'some head';
end
	
	
	
	
i=15;
while i <= 180
	airpar.azimuth = i;
	airpar.channel = 1;
	[h_left,air_info] = load_air(airpar);
	airpar.channel = 0;
	[h_right,air_info] = load_air(airpar);
	irs.left(:,i/15) = h_left;
	irs.right(:,i/15) = h_right;
	irs.room = air_info.room;
    irs.head_azimuth(i/15) = correct_azimuth(rad(360-i+90));
 	irs.apparent_azimuth(i/15) = correct_azimuth(rad(360-i+90));
	irs.apparent_elevation(i/15) = 0;
	i=i+15;
	if i ==30 & irs.distance ~= 3;
		i = i+ 15;
		a = 0
	end
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
outfile = sprintf('%s/AACHEN_%s.mat',outdir,irsset);
save('-v6',outfile,'irs'); 									%v6 is used for the Windows-Version of Octave, better use v7

