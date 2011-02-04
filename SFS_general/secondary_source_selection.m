function ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src)
%SECONDARY_SOURCE_SELECTION selects which secondary sources are active 
%
%   Usage: ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src)
%
%   Input options:
%       x0,y0       - positions of the secondary sources (m)
%       phi         - directions of the secondary sources (rad)
%       xs,ys       - position/direction of the desired source model
%       src         - source type of the virtual source
%                       'pw' - plane wave (xs, ys are the direction of the
%                              plane wave in this case)
%                       'ps' - point source
%                       'fs' - focused source
%       conf        - configuration struct
%
%   Output options:
%       ls_activity - index of the active secondary sources
%
%   SECONDARY_SOURCES_SELECTION(x0,y0,phi,xs,ys,src) returns an vector that
%   indicates which secondary sources are active for the given geometry.
%
%   References:
%       S. Spors, R. Rabenstein, J. Ahrens: "The Theory of Wave Field Synthesis
%       Revisited", in 124th AES Convention, Amsterdam, 2008
%
% see also: secondary_source_positions
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
isargvector({x0,y0,phi},{'x0','y0','phi'});
isargscalar({xs,ys},{'xs','ys'});
isargchar({src},{'src'});


%% ===== Calculation ====================================================

% Direction of secondary sources
nx0 = -sin(phi);
ny0 = cos(phi);

if strcmp('pw',src)
    % === Plane wave ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <n_pw,n_x0> > 0
    % a = <
    %      \ 0, else
    %
    % Direction of plane wave
    nxs = xs / sqrt(xs^2+ys^2);
    nys = ys / sqrt(xs^2+ys^2);
    ls_activity = (( nxs.*nx0 + nys.*ny0 > 0 ));

elseif strcmp('ps',src)
    % === Point source ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <x0-xs,n_x0> > 0
    % a = <
    %      \ 0, else
    %
    ls_activity = (( (x0-xs).*nx0 + (y0-ys).*ny0 > 0 ));

elseif strcmp('fs',src)
    % === Focused source ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <xs-x0,n_x0> > 0
    % a = <
    %      \ 0, else
    ls_activity = (( (xs-x0).*nx0 + (ys-y0).*ny0 > 0 ));
else
    error('%s: %s is not a supported source type!',upper(mfilename),src);
end
