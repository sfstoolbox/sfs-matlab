function kemar2mat(irsset,irspath)
%KEMAR2MAT converts IRs data given by our KEMAR to our mat-file based format
%
%   Usage: kemar2mat(irsset,irspath);
%
%   Input options:
%       irsset  - IR sets measured with VariSphear and KEMAR. Currently the 
%                 following are available:
%                   'RAR_05m'          - HRIR of RAR with 0.5m distance
%                   'RAR_1m'           - HRIR of RAR with 1m distance
%                   'RAR_2m'           - HRIR of RAR with 2m distance
%                   'RAR_3m'           - HRIR of RAR with 3m distance
%                   'RAR3m_small'     - HRIR of RAR with 3m distance and small
%                                         ears
%                   'RAR3m_small2'    - HRIR of RAR with 3m distance and small
%                                         ears. Second reference measurement.
%                   'RAR3m_head_rot'  - HRIR of RAR with 3m distance and only
%                                         head rotation of KEMAR
%                 NOTE: you still have to give the matching path to the given
%                 data set!
%       irspath - path to the directory containing the IR data
%
%   KEMAR2MAT(irsset,irspath) converts the IRs data given by the irsset
%   and stored at the given irspath in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargchar({irsset},{'irsset'});
isargdir({irspath},{'irspath'});


%% ===== Computation =====================================================
% Dir to save the IRs
outdir = 'ir_databases';

% Initialize a new IR struct
irs = new_irs();
% Add common struct entries
irs.fs = 44100;
irs.loudspeaker = 'Genelec 8030A';
irs.room = 'Anechoic chamber ITA TU Berlin';
irs.head_elevation = NaN;
irs.torso_elevation = NaN;
irs.head_position = [0 0 0];
irs.head_reference = [0 1 0];
irs.source_reference = [0 0 0];

if strcmp(irsset,'RAR_05m')
    % irs struct entries
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: torso. Ears: large.'];
    irs.head = 'KEMAR, large ears';
    irs.source_position = [0 0.5 0];
    irs.head_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_0.5m_large_ears';
elseif strcmp(irsset,'RAR_1m')
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: torso. Ears. large.'];
    irs.head = 'KEMAR, large ears';
    irs.source_position = [0 1 0];
    irs.head_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_1m_large_ear';
elseif strcmp(irsset,'RAR_2m')
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: torso. Ears. large.'];
    irs.head = 'KEMAR, large ears';
    irs.source_position = [0 2 0];
    irs.head_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_2m_large_ears_';
elseif strcmp(irsset,'RAR_3m')
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: torso. Ears: large.'];
    irs.head = 'KEMAR, large ears';
    irs.source_position = [0 3 0];
    irs.head_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_3m_large_ears_';
elseif strcmp(irsset,'RAR_3m_small')
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: torso. Ears: small.'];
    irs.head = 'KEMAR, small ears';
    irs.source_position = [0 3 0];
    irs.head_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_3m_';
elseif strcmp(irsset,'RAR_3m_small2')
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: torso. Ears: small. ', ...
         'Reference measurement.'];
    irs.head = 'KEMAR, small ears';
    irs.source_position = [0 3 0];
    irs.head_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_3m_control_';
elseif strcmp(irsset,'RAR_3m_head_rot')
    irs.description = ...
        ['KEMAR measurement in RAR of the TU Berlin. Used elevation ', ...
         'angle: 0°; azimuth resolution: 1°. Rotation: head. Ears: large'];
    irs.head = 'KEMAR, large ears';
    irs.source_position = [0 3 0];
    irs.torso_azimuth = NaN;
    irfilebase = 'KEMAR_1deg_3m_large_ears_head_rotation_';
else
    error('%s: the given irsset is not available.',upper(mfilename));
end

% Calculate the distance between head and source
irs.distance = norm(irs.head_position-irs.source_position);

% Read the data
for ii = 1:360

    if strcmp(irsset,'RAR_3m_head_rot')
        irfile = sprintf('%s/%s%03.0f_%i.mat',irspath,irfilebase,ii,ii-181);
        irs.head_azimuth(ii) = correct_azimuth(-(180-ii+1)/180*pi);
    else
        irfile = sprintf('%s/%s%03.0f_%i.mat',irspath,irfilebase,ii,ii-1);
        irs.torso_azimuth(ii) = correct_azimuth(-(180-ii+1)/180*pi);
    end

    load(irfile);
    irs.apparent_azimuth(ii) = correct_azimuth((180-ii+1)/180*pi);
    irs.apparent_elevation(ii) = correct_azimuth(0);
    irs.left(:,ii) = vspolardata.ir_ch1;
    irs.right(:,ii) = vspolardata.ir_ch2;

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
outfile = sprintf('%s/KEMAR_%s.mat',outdir,irsset);
save('-v7',outfile,'irs');
