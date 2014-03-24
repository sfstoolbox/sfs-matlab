function [x0,idx] = secondary_source_selection(x0,xs,src)
%SECONDARY_SOURCE_SELECTION selects which secondary sources are active
%
%   Usage: [x0,idx] = secondary_source_selection(x0,xs,src)
%
%   Input options:
%       x0          - secondary source positions, directions and weights / m [nx7]
%       xs          - position and for focused sources also direction of the
%                     desired source model / m [1x3] or [1x6]
%       src         - source type of the virtual source
%                       'pw' - plane wave (xs is the direction of the
%                              plane wave in this case)
%                       'ps' - point source
%                       'fs' - focused source
%
%   Output options:
%       x0          - secondary sources / m, containing only the active
%                     ones [mx7]
%       idx         - index of the selected sources from the original x0
%                     matrix [mx1]
%
%   SECONDARY_SOURCE_SELECTION(x0,xs,src) returns only the active secondary
%   sources for the given geometry and virtual source. In addition the index of
%   the chosen secondary sources is returned.
%
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin
%
% see also: secondary_source_positions, secondary_source_tapering

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


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargchar(src);
if strcmp('fs',src) && size(xs,2)~=6
    error(['%s: you have chosen "fs" as source type, then xs has ', ...
        'to be [1x6] including the direction of the focused source.'], ...
        upper(mfilename));
elseif ~strcmp('fs',src) && size(xs,2)~=3
    error(['%s: for all source types beside "fs", the size of xs ', ...
        'has to be [1x3].'],upper(mfilename));
end


%% ===== Calculation ====================================================
% Position and direction of secondary sources
x0_tmp = x0;
nx0 = x0(:,4:6);
x0 = x0(:,1:3);
% Make the size of the xs position the same as the number of secondary sources
% in order to allow x0-xs
xs = repmat(xs,size(x0,1),1);

if strcmp('pw',src)
    % === Plane wave ===
    % direction of the plane wave
    nk = bsxfun(@rdivide,xs,vector_norm(xs,2));
    % secondary source selection
    %
    %      / 1, if nk nx0 > 0
    % a = <
    %      \ 0, else
    %
    % see Wierstorf (2014), p.25 (2.47)
    %
    % Direction of plane wave (nxs) is set above
    idx = (( vector_product(nk,nx0,2)>=eps ));
    x0 = x0_tmp(idx,:);

elseif strcmp('ps',src) || strcmp('ls',src)
    % === Point source ===
    % secondary source selection
    %
    %      / 1, if (x0-xs) nx0 > 0
    % a = <
    %      \ 0, else
    %
    % see Wierstorf (2014), p.26 (2.54) and p.27 (2.59)
    %
    idx = (( vector_product(x0-xs,nx0,2)>=eps ));
    x0 = x0_tmp(idx,:);

elseif strcmp('fs',src)
    % === Focused source ===
    % secondary source selection
    % NOTE: (xs-x0) nx0 > 0 is always true for a focused source
    %
    %      / 1, nxs (xs-x0) > 0
    % a = <
    %      \ 0, else
    %
    % see Wierstorf (2014), p.27 (2.67)
    nxs = xs(:,4:6);
    xs = xs(:,1:3);
    idx = (( vector_product(nxs,xs-x0,2)>=eps ));
    x0 = x0_tmp(idx,:);
else
    error('%s: %s is not a supported source type!',upper(mfilename),src);
end

if size(x0,1)==0
    warning('%s: 0 secondary sources were selected.',upper(mfilename));
end
