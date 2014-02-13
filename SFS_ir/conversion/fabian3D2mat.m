function fabian3D2mat(irsset)
%FABIAN3D2mat converts 3D IRs data given by the FABIAN HRIR measurement  
% from Oldenburg to our mat-file based format
%
%   Usage: fabian3D2mat(irsset) with irsset = load('az0.mat');
%
%   Note: irsset has to be an *.mat file containing an irs-struct 
% 
%   Input options:
%       irsset  - IR sets measured with FABIAN. Currently the 
%                 following are available:
%					'az0.mat'  
%
%   FABIAN3D2mat(irsset) converts the IRs data given by the irsset
%   in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format

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
outdir = 'ir_databases';


%% ===== Computation =====================================================
irs = irsset;
irs.apparent_azimuth = correct_azimuth(rad(irs.source_azimuth)');
irs.apparent_elevation = correct_elevation(rad(irs.source_elevation)');
irs.distance = irs.distance';
irs = rmfield(irs,'source_azimuth');
irs = rmfield(irs,'source_elevation');
irs.ears = 'ear';
irs.head_position = [0;0;0];
irs.head_reference = [0;1;0];

irs.source_reference = [0;0;0];
irs.head_elevation = NaN;
irs.torso_azimuth = NaN;
irs.torso_elevation = NaN;
[x,y,z] =sph2cart(irs.apparent_azimuth,irs.apparent_elevation,irs.distance);
irs.source_position = [x;y;z];

irs = rmfield(irs,'microphone');
irs = rmfield(irs,'itd');

irs = order_irs_fields(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
OutputName = 'FABIAN_3D_anechoic';
save([OutputName '.mat'], 'irs');
