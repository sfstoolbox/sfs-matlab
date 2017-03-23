function ir = sofa_get_data_fire(sofa,idxM,idxE)
%SOFA_GET_DATA_FIRE returns impulse responses from a SOFA file or struct
%
%   Usage: ir = sofa_get_data_fire(sofa,[idxM],[idxE])
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%       idxM    - index of the measurements for which the single impulse
%                 responses should be returned.
%                 idxM could be a single value, then only one impulse response
%                 will be returned, or it can be a vector then all impulse
%                 responses for the corresponding index positions will be
%                 returned (default: all measurements).
%       idxE    - index of the emitter for which the single impulse
%                 responses should be returned (default: all measurements).
%
%   Output parameters:
%       ir      - impulse response (M,2,E,N), where
%                   M ... number of impulse responses
%                   E ... number of emitters (loudspeakers)
%                   N ... samples
%
%   SOFA_GET_DATA_FIRE(sofa,idxM,idxE) returns impulse response of the given
%   SOFA file or struct, specified by idxM and idxE, where idxM defines the
%   measurements and idxE the emitters for which impulse responses should be
%   returned. If no index is specified all data contained in sofa is returned.
%   For the struct the SOFA file has to loaded before with SOFAload().
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   see also: sofa_get_data_fir, sofa_get_header, get_ir, SOFAload

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
nargmax = 3;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    idxE = ':';
elseif nargin==nargmax-2
    idxE = ':';
    idxM = ':';
end


%% ===== Computation ====================================================
if sofa_is_file(sofa)
    header = sofa_get_header(sofa);
    if isnumeric(idxE) && isnumeric(idxM)
        ir = zeros(length(idxM),header.API.R,length(idxE),header.API.N);
        for ii=1:length(idxM)
            for jj=1:length(idxE)
                tmp = SOFAload(sofa,[idxM(ii) 1],'M',[idxE(jj) 1],'E');
                ir(ii,:,jj,:) = tmp.Data.IR;
            end
        end
    elseif isnumeric(idxE)
        ir = zeros(header.API.M,header.API.R,length(idxE),header.API.N);
        for jj=1:length(idxE)
            tmp = SOFAload(sofa,[idxE(jj) 1],'E');
            ir(:,:,jj,:) = tmp.Data.IR;
        end
    elseif isnumeric(idxM)
        ir = zeros(length(idxM),header.API.R,header.API.E,header.API.N);
        for ii=1:length(idxM)
            tmp = SOFAload(sofa,[idxM(ii) 1],'M');
            ir(ii,:,:,:) = tmp.Data.IR;
        end
    else
        tmp = SOFAload(sofa);
        ir = tmp.Data.IR;
    end
else
    ir = sofa.Data.IR(idxM,:,idxE,:);
end
