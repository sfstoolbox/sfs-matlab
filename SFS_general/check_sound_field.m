function check_sound_field(P,t)
%CHECK_SOUND_FIELD checks if we have any activity in the sound field and returns a
%   warning otherwise.
%
%   Usage: check_sound_field(P,t)
%
%   Input parameters:
%       P       - sound field
%       t       - time t / samples
%
%   CHECK_SOUND_FIELD(P,t) checks if the sound field is different from zero.
%   If this is not the case it returns a warning.
%
%   see also: sound_field_imp_wfs, norm_sound_field

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(P);
isargscalar(t);


%% ===== Computation =====================================================
if length(P)==1, return; end
%
if max(abs(P(:)))==0 || all(isnan(P(:)))
    warning('SFS:check_sound_field',...
        ['The activity in the simulated sound field is zero. ',...
         'Maybe you should use another time frame t than %i. '],t);
end
