function fabian2mat(irsset,irspath)
%FABIAN2MAT converts IRs data given by FABIAN to our mat-file based format
%
%   Usage: fabian2mat(irsset,irspath);
%
%   Input options:
%       irsset  - IR sets measured with FABIAN. Currently the following are
%                 available:
%                   'RAR'         - HRIR of RAr of TU Berlin
%                   'Sputnik1'    - BRIR of Sputnik with source at 90°
%                   'Sputnik2'    - BRIR of Sputnik with source at 47°
%                   'Sputnik3'    - BRIR of Sputnik with source at 24°
%                   'Sputnik4'    - BRIR of Sputnik with source at 11°
%                   'Sputnik5'    - BRIR of Sputnik with source at -12°
%                   'Sputnik6'    - BRIR of Sputnik with source at -33°
%                   'Sputnik7'    - BRIR of Sputnik with source at -56°
%                   'Sputnik8'    - BRIR of Sputnik with source at -90°
%                   'audimax_front'
%                   'audimax_back'
%                   'burgtheater_front'
%                   'burgtheater_back'
%                   'friedrichstadtpalast_front'
%                   'friedrichstadtpalast_back'
%                   'udk_kammersaal_front'
%                   'udk_kammersaal_back'
%                   'studio_sweet_spot'
%                 NOTE: you still have to give the matching path to the given
%                 data set!
%       irspath - path to the directory containing the IR data
%
%   FABIAN2MAT(irsset,irspath) converts the IRs data given by the irsset
%   and stored at the given irspath in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format

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
% Common fields
irs.head = 'FABIAN';
irs.ears = 'FABIAN';
irs.fs = 44100;
irs.head_position = [0 0 0];
irs.head_reference = [0 1 0];
irs.source_reference = [0 0 0];
irs.head_elevation = NaN;
irs.torso_elevation = NaN;

if strcmp(irsset,'RAR')
    irs.description = ...
        ['Alex Lindaus anechoic measurements with FABIAN in RAR of the ',...
         'TU berlin. Used elevation angle: 0°; azimuth resolution: ', ...
         '1.25°. NOTE: the real resolution was 2.5°; 1.25° has been ', ...
         'realized by turning the dummy head one time per hand. Rotation: ',...
         'torso.'];
    irs.room = 'Anechoic chamber of ITA TU Berlin';
    irs.source = 'Genelec ?';
    irs.distance = 2.5;
    irs.source_position = [0 2.5 0];
    irs.head_azimuth = NaN;
    angle1 = -180;
    angle2 = 180-1.25;
    step = 1.25;
    irlength = 11025;
elseif strcmp(irsset,'Sputnik1')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 90 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 1.03;
    irs.source_position = [-1.03 0 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik2')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 47 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 1.96;
    irs.source_position = [-1.4335 1.3367 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik3')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 24 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 2.27;
    irs.source_position = [-0.92329 2.0737 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik4')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 11 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 2.86;
    irs.source_position = [-0.54571 2.8075 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik5')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at -12 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 2.92;
    irs.source_position = [0.60710 2.8562 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik6')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at -33 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 2.41;
    irs.source_position = [1.3126 2.0212 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik7')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at -56 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 2.02;
    irs.source_position = [1.6747 1.1296 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik8')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at -90 deg'];
    irs.room = 'Sputnik of T-Labs Berlin';
    irs.source = 'Fostex';
    irs.distance = 1.04;
    irs.source_position = [1.04 0 0];
    irs.torso_azimuth = NaN;
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'audimax_front')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Audimax. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Audimax of T-Labs Berlin';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -80;
    angle2 = 80;
    step = 1;
    irlength = 262144;
elseif strcmp(irsset,'audimax_back')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Audimax. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Audimax of T-Labs Berlin';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -80;
    angle2 = 80;
    step = 1;
    irlength = 262144;
elseif strcmp(irsset,'burgtheater_front')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Burgtheater. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Burgtheater';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -76;
    angle2 = 76;
    step = 2;
    irlength = 262144;
elseif strcmp(irsset,'burgtheater_back')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Burgtheater. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Burgtheater';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -76;
    angle2 = 76;
    step = 2;
    irlength = 262144;
elseif strcmp(irsset,'friedrichstadtpalast_front')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Friedrichstadtpalast. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Friedrichstadtpalast';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -76;
    angle2 = 76;
    step = 2;
    irlength = 262144;
elseif strcmp(irsset,'friedrichstadtpalast_back')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Friedrichstadtpalast. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Friedrichstadtpalast';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -76;
    angle2 = 76;
    step = 2;
    irlength = 262144;
elseif strcmp(irsset,'udk_kammersaal_front')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Kammersaal of UdK. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'UdK Kammersaal';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -76;
    angle2 = 76;
    step = 2;
    irlength = 65536;
elseif strcmp(irsset,'udk_kammersaal_back')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Kammersaal of UdK. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'UdK Kammersaal';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -76;
    angle2 = 76;
    step = 2;
    irlength = 65536;
elseif strcmp(irsset,'studio_sweet_spot')
    irs.description = ...
        ['Alex Lindaus measurements with FABIAN in Studio of TU Berlin. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head. Source at 0 deg'];
    irs.room = 'Studio';
    irs.source = 'Meyer UPL-1';
    irs.distance = 1;
    irs.source_position = [0 1 0];
    irs.torso_azimuth = NaN;
    angle1 = -75;
    angle2 = 75;
    step = 1;
    irlength = 33450;
else
    error('%s: the given irsset is not available.',upper(mfilename));
end

% Read the data
idx = 1;
for phi = angle1:step:angle2

    if strcmp(irsset,'RAR')
        if phi<0
            irfile = sprintf('%s/N%05.0f.wav',irspath,-100*phi);
        else
            irfile = sprintf('%s/P%05.0f.wav',irspath,100*phi);
        end
        irs.torso_azimuth(idx) = -correct_azimuth(rad(phi));
    else
        if phi<0
            irfile = sprintf('%s/P00_N%03.0f.wav',irspath,-10*phi);
        else
            irfile = sprintf('%s/P00_P%03.0f.wav',irspath,10*phi);
        end
        irs.head_azimuth(idx) = -correct_azimuth(rad(phi));
    end

    ir = wavread(irfile);

    irs.apparent_elevation(idx) = correct_elevation(0);
    irs.apparent_azimuth(idx) = ...
        correct_azimuth(rad(phi)-atan2(irs.source_position(1),...
        irs.source_position(2)));
    irs.left(:,idx) = ir(1:irlength,1);
    irs.right(:,idx) = ir(1:irlength,2);

    idx = idx+1;

end

check_irs(irs);
% Reorder entries
irs = correct_irs_angle_order(irs);
irs = order_irs_fields(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
outfile = sprintf('%s/TU_FABIAN_%s.mat',outdir,irsset);
save('-v7',outfile,'irs');
