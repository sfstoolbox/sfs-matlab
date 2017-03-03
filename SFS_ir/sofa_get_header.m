function header = sofa_get_header(sofa)
%SOFA_GET_HEADER returns the header of a SOFA file or struct
%
%   Usage: header = sofa_get_header(sofa)
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%
%   Output parameters:
%       header  - SOFA header
%
%   SOFA_GET_HEADER(sofa) returns the header of the given SOFA file or struct.
%   For the struct the SOFA file has to loaded before with SOFAload().
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   See also: sofa_get_data, sofa_is_file, get_ir, SOFAload 

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


%% ===== Computation ====================================================
% Get only the metadata of the SOFA data set
if sofa_is_file(sofa)
    header = SOFAload(sofa,'nodata');
else
    header = sofa;
    if isfield(sofa.Data,'IR')
        header.Data = rmfield(sofa.Data,'IR');
    end
end
