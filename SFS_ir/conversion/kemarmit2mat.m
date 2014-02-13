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
%   https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
%
%   The MIT HRTF data base can be found here:
%   http://sound.media.mit.edu/resources/KEMAR.html

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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
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
