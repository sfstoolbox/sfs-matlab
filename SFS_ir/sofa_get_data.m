function ir = sofa_get_data(sofa,idx)
%SOFA_GET_DATA returns a single impulse response from a SOFA file or struct
%
%   Usage: ir = sofa_get_data(sofa,idx)
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%       idx     - index of the single impulse response that should be returned
%
%   Output parameters:
%       ir      - impulse response (nx2)
%
%   SOFA_GET_DATA(sofa,idx) returns the single impulse response of the given
%   SOFA file or struct, specified by idx.
%   For the struct the SOFA file has to loaded before with SOFAload().
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   see also: sofa_get_header, get_ir, SOFAload

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


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax)


%% ===== Computation ====================================================
% Check if we have a file or struct
if ~isstruct(sofa) && exist(sofa,'file')
    file = sofa;
    % FIXME: this is currently not working under Octave, because it will
    % always return [1 idx] data
    sofa = SOFAload(file,[idx 1]);
    ir = squeeze(sofa.Data.IR)';
elseif isstruct(sofa) && isfield(sofa.Data,'IR')
    ir = sofa.Data.IR(idx,:,:);
    ir = squeeze(ir)';
else
    error('%s: sofa has to be a file or a SOFA struct.',upper(mfilename));
end
