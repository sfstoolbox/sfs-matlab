function sofa = dummy_irs(nsamples,conf)
% DUMMY_IRS creates a dummy dirac pulse impulse response set
%
%   Usage: irs = dummy_irs([nsamples],[conf])
%
%   Input parameters:
%       nsamples  - length of impulse response in samples, default: 1024
%
%   Output parameters:
%       sofa      - sofa struct
%
%   DUMMY_IRS(nsamples) creates a dummy impulse response data set (Dirac
%   impulse) to check processing without real impulse responses. It has a
%   resolution of 1 deg for phi and theta, its length is given by nsamples.
%
%   See also: SOFAgetConventions

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking input parameters =======================================
nargmin = 0;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin==2
    isargpositivescalar(nsamples);
    isargstruct(conf);
elseif nargin==1
    if isstruct(nsamples)
        conf = nsamples;
        nsamples = 1024;
    else
        conf = SFS_config;
        isargpositivescalar(nsamples);
    end
else
    nsamples = 1024;
    conf = SFS_config;
end


%% ===== Configuration ===================================================
fs = conf.fs;
c = conf.c;
dirac_position = 300;


%% ===== Computation =====================================================
% Angles of dummy irs (in 1deg)
theta = -90:89;
phi = -180:179;
M = length(phi)*length(theta);
ir = zeros(M,2,nsamples);
% Create dirac pulse
ir(:,:,dirac_position) = 1;
% Store data
sofa = SOFAgetConventions('SimpleFreeFieldHRIR');
sofa.Data.IR = ir;
sofa.Data.SamplingRate = fs;
% Add metadata
sofa.GLOBAL_ListenerShortName = 'dummy';
sofa.GLOBAL_History='Created by Sound Field Synthesis Toolbox';
sofa.GLOBAL_Comment = ['HRIR dummy set (Dirac pulse) for testing your',...
                      'frequency response, etc.'];
sofa.ListenerPosition = [0 0 0];
sofa.ListenerView = [1 0 0];
sofa.ListenerUp = [0 0 1];
distance = dirac_position/fs*c;
% Get angles and distance in to the correct format
azimuth = repmat(phi,length(theta),1);
azimuth = azimuth(:);
elevation = repmat(theta,1,length(phi))';
distance = distance.*ones(M,1);
sofa.SourcePosition = [nav2sph(azimuth) elevation distance];
sofa = SOFAupdateDimensions(sofa);
