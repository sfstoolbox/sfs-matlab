function fabianpinta2mat(irsset,irspath)
%FABIANPINTA2MAT converts BRIRs from FABIAN in Pinta to our mat-file based format
%
%   Usage: fabianpinta2mat(irsset,irspath);
%
%   Input options:
%       irsset  - IR sets measured with FABIAN. Currently the following are
%                 available:
%                   'center'      - BRIRs for center position
%                   'center_akg'  - BRIRs for center position with AKG
%                                   headphones mounted on FABIAN
%                   'center_stax' - BRIRs for center position with STAX
%                                   headphones mounted on FABIAN
%                   'front'       - BRIRs for frontal position
%                   'side'        - BRIRs for lateral position
%                 NOTE: you still have to give the matching path to the given
%                 data set!
%       irspath - path to the directory containing the IR data
%
%   FABIANPINTA2MAT(irsset,irspath) converts the IRs data given by the irsset
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
irs.source = 'Elac';
irs.room = 'Pinta at T-Labs Berlin';
irs.fs = 44100;
irs.head_position = [0 0 0];
irs.head_reference = [0 1 0];
irs.source_reference = [0 0 0];
irs.head_elevation = NaN;
irs.torso_elevation = NaN;
irs.torso_azimuth = NaN;
irlength = 180;
angle1 = -80;
angle2 = 80;
step = 1;

if strcmp(irsset,'center')
    irs.description = ...
        ['Alexander Lindaus measurements with FABIAN in Pinta at the ',...
         'T-Labs, central position.', ...
         'Used elevation angle: 0°; azimuth resolution: ', ...
         '1°. Rotation: torso.'];
    itd_offset = -2;
    ild_offset = -2.2989;
elseif strcmp(irsset,'front')
    irs.description = ...
        ['Alexander Lindaus measurements with FABIAN in Pinta at the ',...
         'T-Labs, frontal position.', ...
         'Used elevation angle: 0°; azimuth resolution: ', ...
         '1°. Rotation: torso.'];
    irs.head_position = [0 0.745 0];
    itd_offset = -1;
    ild_offset = -1.3448;
elseif strcmp(irsset,'side')
    irs.description = ...
        ['Alexander Lindaus measurements with FABIAN in Pinta at the ',...
         'T-Labs, side position.', ...
         'Used elevation angle: 0°; azimuth resolution: ', ...
         '1°. Rotation: torso.'];
    irs.head_position = [0.755 0 0];
    itd_offset = -1;
    ild_offset = -1.6388;
elseif strcmp(irsset,'center_akg')
    irs.description = ...
        ['Alexander Lindaus measurements with FABIAN in Pinta at the ',...
         'T-Labs, central position with AKG headphones.', ...
         'Used elevation angle: 0°; azimuth resolution: ', ...
         '1°. Rotation: torso.'];
    irs.ears = 'FABIAN + AKG K 601 headphones';
    itd_offset = 0;
    ild_offset = -1.1820;
elseif strcmp(irsset,'center_stax')
    irs.description = ...
        ['Alexander Lindaus measurements with FABIAN in Pinta at the ',...
         'T-Labs, central position with STAX headphones.', ...
         'Used elevation angle: 0°; azimuth resolution: ', ...
         '1°. Rotation: torso.'];
    irs.ears = 'FABIAN + STAX headphones';
    itd_offset = 0;
    ild_offset = -1.5405;
else
    error('%s: the given irsset is not available.',upper(mfilename));
end

angles = rad(0:360/56:355);
for ii = 1:56
    irs.source_position = [1.5*sin(angles(ii)) 1.5*cos(angles(ii)) 0];
    irs.distance = norm(irs.head_position-irs.source_position);
    name = sprintf('%s_src%i',irsset,ii);
    path = sprintf('%s/source%i',irspath,ii);
    % Read the data
    idx = 1;
    for phi = angle1:step:angle2

        if phi<0
            irfile = sprintf('%s/P00_N%03.0f.wav',path,-10*phi);
        else
            irfile = sprintf('%s/P00_P%03.0f.wav',path,10*phi);
        end
        irs.head_azimuth(idx) = correct_azimuth(-rad(phi+itd_offset));

        ir = wavread(irfile);

        irs.apparent_elevation(idx) = correct_elevation(0);
        irs.apparent_azimuth(idx) = -correct_azimuth(-rad(phi+itd_offset)+...
            atan2(irs.source_position(1)-irs.head_position(1),...
                  irs.source_position(2)-irs.head_position(2)));
        irs.left(:,idx) = gaindb(ir(1:irlength,1),ild_offset);
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
    outfile = sprintf('%s/TU_FABIAN_pinta_%s.mat',outdir,name);
    save('-v7',outfile,'irs');

end

% Create different sets by using a fixed head position and the different sources
% as "simulated" head movements.
% NOTE: this kind of data set are of course needed in order to simulate the
% circular loudspeaker array!
irs.head_azimuth = [];
irs.apparent_azimuth = [];
irs.apparent_elevation = [];
irs.left = [];
irs.right = [];
angles = rad(0:360/56:355);
for phi = angle1+abs(itd_offset):step:angle2-abs(itd_offset)
    irs.source_position = [1.5*sin(rad(phi)) 1.5*cos(rad(phi)) 0];
    irs.distance = norm(irs.head_position-irs.source_position);
    name = sprintf('%s_%ideg',irsset,phi);
    idx = 1;
    for ii = 1:56
        source_position = [1.5*-sin(angles(ii)) 1.5*cos(angles(ii)) 0];
        % Read the data
        path = sprintf('%s/source%i',irspath,ii);
        if phi-itd_offset<0
            irfile = sprintf('%s/P00_N%03.0f.wav',path,-10*(phi-itd_offset));
        else
            irfile = sprintf('%s/P00_P%03.0f.wav',path,10*(phi-itd_offset));
        end
        ir = wavread(irfile);
        irs.head_azimuth(idx) = -correct_azimuth(rad(phi));
        irs.apparent_elevation(idx) = correct_elevation(0);
        irs.apparent_azimuth(idx) = correct_azimuth(rad(phi) + ...
            atan2(-(source_position(1)-irs.head_position(1)),...
                  source_position(2)-irs.head_position(2)));
        irs.left(:,idx) = gaindb(ir(1:irlength,1),ild_offset);
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
    outfile = sprintf('%s/TU_FABIAN_pinta_%s.mat',outdir,name);
    save('-v7',outfile,'irs');

end
