function x0 = secondary_source_selection(x0,xs,src,xref)
%SECONDARY_SOURCE_SELECTION selects which secondary sources are active
%
%   Usage: ls_activity = secondary_source_selection(x0,xs,src,[xref])
%
%   Input options:
%       x0          - secondary source positions and directions (m)
%       xs          - position of the desired source model (m)
%       src         - source type of the virtual source
%                       'pw' - plane wave (xs is the direction of the
%                              plane wave in this case)
%                       'ps' - point source
%                       'fs' - focused source
%       xref        - position of reference point (conf.xref). This is needed
%                     for focused sources, and will define the direction of the
%                     focused source
%
%   Output options:
%       x0          - secondary sources, containing only the active ones
%
%   SECONDARY_SOURCES_SELECTION(x0,xs,src) returns only the active secondary
%   sources for the given geometry and virtual source.
%
%   References:
%       S. Spors, R. Rabenstein, J. Ahrens: "The Theory of Wave Field Synthesis
%       Revisited", in 124th AES Convention, Amsterdam, 2008
%
% see also: secondary_source_positions, secondary_source_number, tapering_window

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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
xs = position_vector(xs);
isargchar(src);
if nargin==nargmax
    xref = position_vector(xref);
elseif strcmp('fs',src)
    error(['%s: you have chosen "fs" as source type, then xref is ', ...
        'needed as fourth argument.'],upper(mfilename));
else
    xref = [0,0,0];
end


%% ===== Calculation ====================================================
% Position and direction of secondary sources
x0_tmp = x0;
nx0 = x0(:,4:6);
x0 = x0(:,1:3);
% direction of the plane wave
nxs = xs / norm(xs);
% Make the size of the xs position the same as the number of secondary sources
% in order to allow x0-xs
xs = repmat(xs,size(x0,1),1);
xref = repmat(xref,size(x0,1),1);
nxs = repmat(nxs,size(x0,1),1);

if strcmp('pw',src)
    % === Plane wave ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <n_pw,n_x0> > 0
    % a = <
    %      \ 0, else
    %
    % Direction of plane wave (nxs) is set above
    x0 = x0_tmp(diag(nxs*nx0')>=eps,:);

elseif strcmp('ps',src) || strcmp('ls',src)
    % === Point source ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <x0-xs,n_x0> > 0
    % a = <
    %      \ 0, else
    %
    x0 = x0_tmp(diag((x0-xs)*nx0')>0,:);

elseif strcmp('fs',src)
    % === Focused source ===
    % secondary source selection (Spors 2008)
    % NOTE: <xs-x0,nx0> > 0 is always true for a focused source
    %
    %      / 1, <xs-xref,x_0-xs> > 0
    % a = <
    %      \ 0, else
    %
    x0 = x0_tmp(diag((xs-xref)*(x0-xs)')>0,:);
else
    error('%s: %s is not a supported source type!',upper(mfilename),src);
end
