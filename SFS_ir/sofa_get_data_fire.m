function ir = sofa_get_data_fire(sofa,idxM,idxE)
%SOFA_GET_DATA returns impulse responses from a SOFA file or struct
%
%   Usage: ir = sofa_get_data(sofa,[idx])
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%       idx     - index of the single impulse responses that should be returned
%                 idx could be a single value, then only one impulse response
%                 will be returned, or it can be a vector then all impulse
%                 responses for the corresponding index positions will be
%                 returned.
%                 If no index is specified all data will be returned.
%
%   Output parameters:
%       ir      - impulse response (M,2,E,N), where
%                   M ... number of impulse responses
%                   E ... number of emitters (loudspeakers)
%                   N ... samples
%
%   SOFA_GET_DATA_FIRE(sofa,idxM,idxE) returns impulse response of the given
%   SOFA file or struct, specified by idx. If no idx is specified all data
%   contained in sofa is returned.
%   For the struct the SOFA file has to loaded before with SOFAload().
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   see also: sofa_get_header, sofa_is_file, get_ir, SOFAload

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
