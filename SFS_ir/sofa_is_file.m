function boolean = sofa_is_file(sofa)
%SOFA_IS_FILE returns 1 for a sofa file, 0 for a sofa struct or an error otherwise
%
%   Usage: number = sofa_check(sofa)
%
%   Input parameters:
%       sofa    - sofa struct or file name
%
%   Output parameters:
%       number  - 1: sofa is a file
%                 0: sofa is a struct
%
%   SOFA_CHECK(sofa) checks if the given sofa is a file or a struct. In the
%   first case a 1 is returned, in the second case a 0. If none of both is true
%   an error will be thrown.
%
%   See also: sofa_get_header, sofa_get_data

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
% No checking for speedup
%nargmin = 1;
%nargmax = 1;
%narginchk(nargmin,nargmax)


%% ===== Main ===========================================================
if isstruct(sofa) && isfield(sofa,'GLOBAL_Conventions') && ...
        strcmp('SOFA',sofa.GLOBAL_Conventions)
    boolean = false;
elseif ~isstruct(sofa) && exist(sofa,'file')
    boolean = true;
else
    error('%s: sofa has to be a file or a SOFA struct.',upper(mfilename));
end
