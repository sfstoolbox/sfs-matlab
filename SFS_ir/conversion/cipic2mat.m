function cipic2mat(irsset,irspath)
%CIPIC2MAT converts IRs data given from the CIPIS measurements to our mat-file
%based format
%
%   Usage: cipic2mat(irsset,irspath);
%
%   Input options:
%       irsset  - IR sets measured from CIPIS.
%                 Currently the following are available:
%                   'large_pinna_frontal'  - HRIR from KEMAR with large ears,
%                                            only from the front (spacing of
%                                            2.8125°)
%                   'small_pinna_frontal'  - HRIR from KEMAR with small ears,
%                                            only from the front (spacing of
%                                            2.8125°)
%                   'large_pinna_final'    - HRIR from KEMAR with large ears
%                                            (spacing of 5°)
%                   'small_pinna_final'    - HRIR from KEMAR with small ears
%                                            (spacing of 5°)
%                   'subject_xxx'          - for example subject_003 , here
%                                            real persons were taken, also
%                                            elevation-data available
%                   'subject_021'          - KEMAR with large pinnae, also
%                                            elevation-data available
%                   'subject_165'          - KEMAR with small pinnae, also
%                                            elevation-data available
%                 NOTE: you still have to give the matching path to the
%                 CIPIC-data
%       irspath - path to the directory containing the CIPIC-data
%
%   CIPIC2MAT(irsset,irspath) converts the IRs data given by the irsset
%   and stored at the given irspath in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
%
%   The CIPIC data set can be found here:
%   http://interface.cipic.ucdavis.edu/sound/hrtf.html

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargchar(irsset);
isargdir(irspath);


%% ===== Computation =====================================================
% Dir to save the IRs
outdir = 'ir_databases';

% Initialize a new IR struct
irs = new_irs();

% Add common struct entries
irs.fs = 44100;

irs.head_azimuth = NaN;
irs.torso_azimuth = NaN;
irs.head_position = [0,0,0]';
irs.head_reference = [0 1 0]';
irs.source_reference = [0 0 0]';
irs.head_elevation = NaN;
irs.torso_elevation = NaN;
irs.room = 'CIPIC';
irs.source ='5x BOSE ACOUSTIMASS 6,4cm';
irs.room_corners = [0 0 -1.7;]';


if strcmp(irsset,'large_pinna_frontal')
    irs.head = 'KEMAR';
    irs.description = 'CIPIC measurement with KEMAR,spacing of 2.8125°';
    irs.ears = 'large';
elseif strcmp(irsset,'small_pinna_frontal')
    irs.head = 'KEMAR';
    irs.description = 'CIPIC measurement with KEMAR,spacing of 2.8125°';
    irs.ears = 'small';
elseif strcmp(irsset,'large_pinna_final')
    irs.head = 'KEMAR';
    irs.description = 'CIPIC measurement with KEMAR,spacing of 5°';
    irs.ears = 'large';
elseif strcmp(irsset,'small_pinna_final')
    irs.head = 'KEMAR';
    irs.description = 'CIPIC measurement with KEMAR,spacing of 5°';
    irs.ears = 'small';
elseif strncmp(irsset,'subject',7)
    irs.head = sprintf('%s for more details see anthro.mat',irsset);
    irs.ears = sprintf('%s for more details see anthro.mat',irsset);
    irs.description = ...
        ['CIPIC; azimuth resolution: 5° (-45° to +45°) and +- 55° 65° 80° '];

end



if strcmp(irsset,'large_pinna_frontal') | strcmp(irsset,'small_pinna_frontal')
    irfile = sprintf('%s/special_kemar_hrir/kemar_frontal/%s.mat',irspath,irsset);
    load(irfile);
    for ii = 1:99
        theta = 0;                      % elevation
        phi = rad(135-(ii-1)*2.8125);   % azimuth
        r = 1;                          % radius

        irs.distance(ii) = r;
        [x y z] = sph2cart(phi,theta,r);
        irs.source_position(:,ii) = [x y z]';

        irs.apparent_azimuth(ii) = correct_azimuth(-phi);
        irs.apparent_elevation(ii) = theta;

        irs.left(:,ii) = left(:,ii);
        irs.right(:,ii) = right(:,ii);
    end
elseif strcmp(irsset,'large_pinna_final') | strcmp(irsset,'small_pinna_final')
    irfile = sprintf('%s/special_kemar_hrir/kemar_horizontal/%s.mat',irspath,irsset);
    load(irfile);
    for ii = 1:72
        theta = 0;                  % elevation
        phi = rad(180-(ii-1)*5);    % azimuth
        r = 1;                      % radius

        irs.distance(ii) = r;
        [x y z] = sph2cart(phi,theta,r);
        irs.source_position(:,ii) = [x y z]';

        irs.apparent_azimuth(ii) = correct_azimuth(-phi);
        irs.apparent_elevation(ii) = round(theta*10000)/10000;

        irs.left(:,ii) = left(:,ii);
        irs.right(:,ii) = right(:,ii);
    end

elseif strncmp(irsset,'subject',7)
    irfile = sprintf('%s/standard_hrir_database/%s/hrir_final.mat',irspath,irsset);
    load(irfile);
    azimuths = [-80 -65 -55 -45:5:45 55 65 80];
    elevations = -45:5.625:235;

    ii=1;
    for a = 1:25
        for e = 1:50
            theta = rad(elevations(e));         % elevation
            phi = rad(azimuths(a));             % azimuth
            r = 1;                              % radius


            if theta > pi/2                     % behind should be defined by the azimuth instead of elevation
                theta = pi - theta;
                phi = pi-phi;
            end

            irs.distance(ii) = r;
            [x y z] = sph2cart(phi,theta,r);

            irs.source_position(:,ii) = [x y z]';
            irs.apparent_azimuth(ii) = correct_azimuth(-phi);
            irs.apparent_elevation(ii) = round(theta*100000)/100000;

            irs.left(:,ii) = hrir_l(a,e,:);
            irs.right(:,ii) = hrir_r(a,e,:);
            ii = ii+1;
        end
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
outfile = sprintf('%s/CIPIC_%s.mat',outdir,irsset);
%v6 is used for the Windows-Version of Octave, better use v7
save('-v7',outfile,'irs');
