function ir = sofa_get_data_fir(sofa,idx)
%SOFA_GET_DATA_FIR returns impulse responses from a SOFA file or struct
%
%   Usage: ir = sofa_get_data_fir(sofa,[idx])
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%       idx     - index of the single impulse responses that should be returned.
%                 idx could be a single value, then only one impulse response
%                 will be returned, or it can be a vector then all impulse
%                 responses for the corresponding index positions will be
%                 returned (default: return all impulse responses)
%
%   Output parameters:
%       ir      - impulse response (M,2,N), where
%                   M ... number of impulse responses
%                   N ... samples
%
%   SOFA_GET_DATA_FIR(sofa,idx) returns impulse response of the given
%   SOFA file or struct, specified by idx. If no idx is specified all data
%   contained in sofa is returned.
%   For the struct the SOFA file has to loaded before with SOFAload().
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   See also: sofa_get_data_fire, sofa_get_header, get_ir, SOFAload

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
nargmax = 2;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    idx = [];
else
    isargvector(idx);
end


%% ===== Computation ====================================================
if length(idx)==0
    if sofa_is_file(sofa)
        sofa = SOFAload(sofa);
    end
    ir = sofa.Data.IR;
else
    header = sofa_get_header(sofa);
    if sofa_is_file(sofa)
        ir = zeros(length(idx),2,header.API.N);
        for ii=1:length(idx)
            tmp = SOFAload(sofa,[idx(ii) 1]);
            ir(ii,:,:) = tmp.Data.IR;
        end
    else
        ir = sofa.Data.IR(idx,:,:);
    end
end
