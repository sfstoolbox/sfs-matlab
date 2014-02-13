function bkoldenburg2mat(irsset,irspath)
%BKOLDENBURG2MAT converts IRs from BK dummy head from Oldenburg to irs format
%
%   Usage: bkoldenburg2mat(irspath);
%
%   Input parameters:
%       irsset - IR sets:
%                   'RAR_08m'   - HRIR of RAR with 0.8 m distance.
%                   'RAR_3m'    - HRIR of RAR with 3 m distance.
%       irspath - path to the directory containing the IR data
%
%   BKOLDENBURG2MAT(irsset,irspath) converts the IRs data given by the irsset
%   and measured with a B&K dummy head in Oldenburg and stored at the given
%   irspath in our own mat-file based format. For format details, see:
%   https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
%
%   The Oldenburg HRTF data base can be found here:
%   http://medi.uni-oldenburg.de/hrir/

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
