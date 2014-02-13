function ir = compensate_headphone(ir,conf)
%COMPENSATE_HEADPHONE applies a headphone compensation to the IR
%
%   Usage: ir = compensate_headphone(ir,[conf])
%
%   Input parameters:
%       ir      - Impulse response to which the compensation should be applied
%       conf    - optional configuration struct (see SFS_config)
%
%   Output:
%       ir      - Impulse response which is compensated for the given headphone
%
%   COMPENSATE_HEADPHONE(ir,conf) applies a headphone compensation to the
%   given impulse response. Which headphone compensation it should use is
%   mentioned in the conf struct.
%   The compensation is only applied, if the conf.ir.usehcomp value is not false.
%
%   see also: ir_wfs, ir_point_source, ir_generic

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
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(ir);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
usehcomp = conf.ir.usehcomp;


%% ===== Computation =====================================================
if(usehcomp)
    lenir = size(ir,1);
    % Read headphone compensation filter
    hcomp = wavread(conf.ir.hcompfile);
    % Check if the IR has the right length for the filter
    if lenir<length(hcomp)
        warning(['The length of the used IR is shorter than the headphone ', ...
            'compensation filter.']);
    end
    % Apply filter
    ir = convolution(hcomp,ir);
    ir = fix_length(ir,lenir);
end
