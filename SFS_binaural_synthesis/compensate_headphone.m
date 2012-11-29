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
%   The compensation is only applied, if the conf.usehcomp value is not false.
%
%   see also: SFS_config, ir_wfs_25d, ir_point_source

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
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
usehcomp = conf.usehcomp;
hcomplfile = conf.hcomplfile;
hcomprfile = conf.hcomprfile;


%% ===== Computation =====================================================
if(usehcomp)
    % Read headphone compensation filter
    hcompl = wavread(hcomplfile);
    hcompr = wavread(hcomprfile);
    hcomp = [hcompl hcompr];
    % Check if the IR has the right length for the filter
    if length(ir(:,1))<length(hcompl)
        warning(['The length of the used IR is shorter than the headphone ', ...
            'compensation filter.']);
    end
    % Apply filter
    % The following is the original code from Sascha, but it will work only if
    % the length of the IR is sufficient greater than the length of the
    % headphone compensation filter.
    %ir(:,1) = conv(hcomp(:,1),ir(1:end-length(hcomp)+1,1));
    %ir(:,2) = conv(hcomp(:,2),ir(1:end-length(hcomp)+1,2));
    % Therefore we use this one
    tmp1 = conv(hcomp(:,1),ir(:,1));
    tmp2 = conv(hcomp(:,2),ir(:,2));
    len = length(ir(:,1));
    ir(:,1) = tmp1(1:len);
    ir(:,2) = tmp2(1:len);
end
