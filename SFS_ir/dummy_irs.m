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
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
distance = 1;
dirac_position = 1;


%% ===== Computation =====================================================
ir = zeros(1,2,nsamples);
% Create dirac pulse in first sample as the delay corresponding to the distance
% is handled by get_ir() later on
ir(:,:,dirac_position) = 1/(4*pi);
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
sofa.SourcePosition = [nav2sph(azimuth) elevation distance];
sofa = SOFAupdateDimensions(sofa);
