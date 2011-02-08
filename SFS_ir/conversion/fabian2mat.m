function fabian2mat(irsset,irspath)
%FABIAN2MAT converts IRs data given by FABIAN to our mat-file based format
%
%   Usage: fabian2mat(irsset,irspath);
%
%   Input options:
%       irsset  - IR sets measured with FABIAN. Currently the following are
%                 available:
%                   'RAR'       - HRIR of RAr of TU Berlin
%                   'Sputnik1'  - BRIR of Sputnik with source at 90°
%                   'Sputnik2'  - BRIR of Sputnik with source at 47°
%                   'Sputnik3'  - BRIR of Sputnik with source at 24°
%                   'Sputnik4'  - BRIR of Sputnik with source at 11°
%                   'Sputnik5'  - BRIR of Sputnik with source at -12°
%                   'Sputnik6'  - BRIR of Sputnik with source at -33°
%                   'Sputnik7'  - BRIR of Sputnik with source at -56°
%                   'Sputnik8'  - BRIR of Sputnik with source at -90°
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
% Initialize a new IR struct
irs = {};
if strcmp(irsset,'RAR')
    irs.description = ...
        ['Alex Lindaus anechoic measurements with FABIAN in RAR of the ',...
         'TU berlin. Used elevation angle: 0°; azimuth resolution: ', ...
         '1.25°. NOTE: the real resolution was 2.5°; 1.25° has been ', ...
         'realized by turning the dummy head one time per hand. Rotation: ',...
         'torso.'];
    irs.tag = 'HRIR';
    irs.distance = 2.5;
    irs.direction = 0;
    irs.fs = 44100;
    outdir = 'measurements/HRIRs';
    angle1 = -180;
    angle2 = 180-1.25;
    step = 1.25;
    irlength = 11025;
elseif strcmp(irsset,'Sputnik1')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 1.03;
    irs.direction = 90/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik2')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 1.96;
    irs.direction = 47/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik3')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 2.27;
    irs.direction = 24/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik4')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 2.86;
    irs.direction = 11/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik5')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 2.92;
    irs.direction = -12/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik6')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 2.41;
    irs.direction = -33/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik7')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 2.02;
    irs.direction = -56/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
elseif strcmp(irsset,'Sputnik8')
    irs.description = ...
        ['Alex Lindaus echoic measurements with FABIAN in Sputnik. ', ...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg. ', ...
         'Rotation: head.'];
    irs.tag = 'BRIR';
    irs.distance = 1.04;
    irs.direction = -90/180*pi;
    irs.fs = 44100;
    outdir = 'measurements/BRIRs';
    angle1 = -90;
    angle2 = 90;
    step = 1;
    irlength = 44100;
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
    else
        if phi<0
            irfile = sprintf('%s/P00_N%03.0f.wav',irspath,-10*phi);
        else
            irfile = sprintf('%s/P00_P%03.0f.wav',irspath,10*phi);
        end
    end

    ir = wavread(irfile);

    irs.elevation(idx) = correct_angle(0);
    irs.azimuth(idx) = correct_angle(phi/180*pi);
    irs.left(:,idx) = ir(1:irlength,1);
    irs.right(:,idx) = ir(1:irlength,2);

    idx = idx+1;

end

% Reorder entries
irs = correct_irs_angle_order(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
outfile = sprintf('%s/FABIAN_%s.mat',outdir,irsset);
save('-v7',outfile,'irs');
