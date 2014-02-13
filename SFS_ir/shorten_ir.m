function ir = shorten_ir(ir,nsamples)
%SHORTEN_IR shortens a IR
%
%   Usage: ir = shorten_ir(ir,nsamples)
%
%   Input parameters:
%       ir          - IR signal with length x channels
%       nsamples    - length of the target IR
%
%   Output paramteres:
%       ir          - IR signal with nsamples x n
%
%   SHORTEN_HRIR(ir,nsamples) shortens a given IR to the given number of samples
%   nsamples and applying a 5% long hanning window.
%
%   see also: SFS_config, read_irs, intpol_ir, reduce_ir

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
isargpositivescalar(nsamples);


%% ===== Computation ====================================================

% Window IR
win = hann_window(0,ceil(0.05*nsamples),nsamples);

ir = ir(1:nsamples,:) .* repmat(win,[1 size(ir,2)]);
