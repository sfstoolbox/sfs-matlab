function [ir,delay] = wfs_preequalization(ir,conf)
%WFS_PREEQUALIZATION applies a pre-equalization filter for WFS
%
%   Usage: ir = wfs_preequalization(ir,conf)
%
%   Input parameters:
%       ir      - signal to which the pre-equalization filter should be applied
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - signal with applied pre-equalization
%       delay   - additional delay added by pre-equalization
%
%
%   WFS_PREEQUALIZATION(ir,conf) applies the pre-equalization filter for
%   Wave Field Synthesis to the given impulse response.
%
%   See also: wfs_fir_prefilter, wfs_iir_prefilter, driving_function_imp_wfs

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(ir);
isargstruct(conf);


%% ===== Configuration ==================================================
usehpre = conf.wfs.usehpre;
fs = conf.fs;


%% ===== Computation =====================================================
% Check if we should procide
if ~usehpre
    delay = 0;
    return;
end
% Store original length
len_ir = size(ir,1);
% Get the filter
if strcmp('FIR',conf.wfs.hpretype)
    % Get FIR filter
    hpre = wfs_fir_prefilter(conf);
    % Apply filter
    ir = convolution(hpre,ir);
    % Delay in s added by filter
    delay = conf.wfs.hpreFIRorder/2 / fs;
elseif strcmp('IIR',conf.wfs.hpretype)
    if len_ir == 1
        % Happens when called from driving_function_imp_wfs()
        % Zeropadding to length conf.N
       ir = [ir; zeros(conf.N-1,size(ir,2))];
    end
    % Get IIR filter
    hpre = wfs_iir_prefilter(conf);
    % Apply filter
    ir = sosfilt(hpre.sos,ir,1);
    % IIR is minimum phase, so no proper delay introduced
    delay = 0;
else
    error('%s: %s is an unknown filter type.',upper(mfilename),hpretype);
end
% Correct length of ir
if strcmp('FIR',conf.wfs.hpretype) && (len_ir>length(hpre)+1)
    ir = ir(1:len_ir,:);
end
