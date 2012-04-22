function brir = compensate_headphone(brir,conf)
%COMPENSATE_HEADPHONE applies a headphone compensation to the BRIR
%
%   Usage: brir = compensate_headphone(brir,[conf])
%
%   Input parameters:
%       brir    - BRIR to which the compensation should be applied
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       brir    - BRIR which is compensated for the given headphone
%
%   COMPENSATE_HEADPHONE(brir,conf) applies a headphone compensation to the
%   given BRIR. Which headphone compensation it should use is mentioned in the
%   conf struct.
%   The compensation is only applied, if the conf.usehcomp value is not false.
%
%   see also: SFS_config, brs_wfs_25d

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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox      sfs-toolbox@gmail.com *
%*****************************************************************************

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(brir);
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
    % Check if the BRIR has the right length for the filter
    if length(brir(:,1))<length(hcompl)
        warning(['The length of the used IR is shorter than the headphone ', ...
            'compensation filter.']);
    end
    % Apply filter
    % The following is the original code from Sascha, but it will work only if
    % the length of the IR is sufficient greater than the length of the
    % headphone compensation filter.
    %brir(:,1) = conv(hcomp(:,1),brir(1:end-length(hcomp)+1,1));
    %brir(:,2) = conv(hcomp(:,2),brir(1:end-length(hcomp)+1,2));
    % Therefore we use this one
    tmp1 = conv(hcomp(:,1),brir(:,1));
    tmp2 = conv(hcomp(:,2),brir(:,2));
    len = length(brir(:,1));
    brir(:,1) = tmp1(1:len);
    brir(:,2) = tmp2(1:len);
end
