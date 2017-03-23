function boolean = sofa_is_file(sofa)
%SOFA_CHECK returns 1 for a sofa file, 0 for a sofa struct or an error otherwise
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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax)


%% ===== Main ===========================================================
if ~isstruct(sofa) && exist(sofa,'file')
    boolean = true;
elseif isstruct(sofa) && isfield(sofa,'GLOBAL_Conventions') && ...
       strcmp('SOFA',sofa.GLOBAL_Conventions)
    boolean = false;
else
    error('%s: sofa has to be a file or a SOFA struct.',upper(mfilename));
end
