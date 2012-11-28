function check_wave_field(P,frame)
%CHECK_WAVE_FIELD checks if we have any activity in the wave field and returns a
%   warning otherwise.
%
%   Usage: check_wave_field(P,frame)
%
%   Input parameters:
%       P       - wave field
%       frame   - used time frame
%
%   CHECK_WAVE_FIELD(P,frame) checks if the wave field is different from zero.
%   If this is not the case it returns a warning.
%
%   see also: wave_field_imp_wfs_25d, norm_wave_field

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(P);
isargscalar(frame);


%% ===== Computation =====================================================
if max(abs(P(:)))==0 || all(isnan(P(:)))
    warning('SFS:check_wave_field',...
        ['The activity in the simulated wave field is zero. ',...
         'Maybe you should use another time frame than %i. ', ...
         'You can set the time frame with conf.frame.'],frame);
end
