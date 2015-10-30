function P = norm_sound_field(P,conf)
%NORM_SOUND_FIELD normalizes the sound field
%
%   Usage: P = norm_sound_field(P,conf)
%
%   Input options:
%       P       - sound field
%       conf    - configuration struct (see SFS_config)
%
%   Output options:
%       P       - normalized sound field
%
%   NORM_SOUND_FIELD(P,conf) normalizes the given sound field P. This depends on
%   the conf.plot.normalisation setting. It can be one of the following:
%       'auto'   - if the given absolute sound field value at the center is
%                  > 0.3 it uses automatically 'center', otherwise it uses 'max'
%       'center' - normalises to center of sound field == 1
%       'max'    - normalises to max of sound field == 1
%
%   See also: plot_sound_field

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargnumeric(P);
isargstruct(conf);


%% ===== Configuration ===================================================
method = conf.plot.normalisation;


%% ===== Computation =====================================================
if strcmp('auto',method)
    % If sound field at center > 0.3 normalise to center sound field == 1,
    % otherwise normalise to max sound field == 1
    if abs(P(round(end/2),round(end/2)))>0.3
        method = 'center';
    else
        method = 'max';
    end
end
if strcmp('center',method)
    % Center of sound field == 1
    P = P/max(abs(P(round(end/2),round(end/2))));
elseif strcmp('max',method)
    % Max of sound field == 1
    P = P/max(abs(P(:)));
else
    error(['%s: conf.plot.normalisation has to be ''auto'', ''center'' or ', ...
           '''max''.'],upper(mfilename));
end
