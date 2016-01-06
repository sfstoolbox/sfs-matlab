function sofa = dummy_irs(nsamples,conf)
% DUMMY_IRS creates a dummy dirac pulse impulse response set
%
%   Usage: irs = dummy_irs([nsamples],conf)
%
%   Input parameters:
%       nsamples  - length of impulse response in samples, default: 1024
%       conf      - configuration struct (see SFS_config)
%
%   Output parameters:
%       sofa      - sofa struct
%
%   DUMMY_IRS(nsamples,conf) creates a dummy impulse response data set (Dirac
%   impulse) to check processing without real impulse responses. It returns only
%   one Dirac impulse, which is then applied for all direction if you for
%   example use it together with ir_wfs().
%
%   See also: SOFAgetConventions, get_ir, ir_wfs

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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = nsamples;
    nsamples = 1024;
end
isargpositivescalar(nsamples);
isargstruct(conf);


%% ===== Configuration ===================================================
fs = conf.fs;
c = conf.c;
dirac_position = 300;


%% ===== Computation =====================================================
ir = zeros(1,2,nsamples);
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
elevation = 0;
azimuth = 0;
distance = dirac_position/fs*c;
sofa.SourcePosition = [nav2sph(azimuth) elevation distance];
sofa = SOFAupdateDimensions(sofa);
