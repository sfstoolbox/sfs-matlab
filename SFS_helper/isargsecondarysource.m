function isargsecondarysource(varargin)
%ISARGSECONDARYSOURCE tests if the given arg is a matrix containing secondary
%   source positions and directions and return an error otherwise
%
%   Usage: isargsecondarysource(args)
%
%   Input options:
%       args        - list of args
%
%   ISARGSECONDARYSOURCE(args) tests if all given args are a matrix
%   containing secondary source positions and directions and returns
%   an error otherwise.
%
%   See also: isargmatrix

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


%% ===== Checking for vector =============================================
for ii = 1:nargin
   
    x0 = varargin{ii};
    
    if ~isnumeric(x0) || ndims(x0)~=2 || size(x0,2)~=7 
        error(['%s need to be a nx7 matrix containing the ', ...
            'secondary sources positions and directions ', ...
            'and weights.'], ...
            inputname(ii));
    end
    if ~all(abs(vector_norm(x0(:,4:6),2)-1)<1e-10)
        error(['The norm of the direction of %s is not 1. ', ...
            'You can use the direction_vector function to get ', ...
            'correct direction values.'],inputname(ii));
    end
end
