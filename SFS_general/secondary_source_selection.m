function ls_activity = secondary_source_selection(x0,xs,src)
%SECONDARY_SOURCE_SELECTION selects which secondary sources are active
%
%   Usage: ls_activity = secondary_source_selection(x0,xs,src)
%
%   Input options:
%       x0          - secondary source positions and directions (m)
%       xs          - position of the desired source model (m)
%       src         - source type of the virtual source
%                       'pw' - plane wave (xs is the direction of the
%                              plane wave in this case)
%                       'ps' - point source
%                       'fs' - focused source
%       conf        - configuration struct
%
%   Output options:
%       ls_activity - index of the active secondary sources
%
%   SECONDARY_SOURCES_SELECTION(x0,xs,src) returns an vector that
%   indicates which secondary sources are active for the given geometry.
%
%   References:
%       S. Spors, R. Rabenstein, J. Ahrens: "The Theory of Wave Field Synthesis
%       Revisited", in 124th AES Convention, Amsterdam, 2008
%
% see also: secondary_source_positions, secondary_source_number, tapwin
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargsecondarysource(x0)
isargposition(xs);
xs = position_vector(xs);
isargchar(src);


%% ===== Calculation ====================================================

% Position and direction of secondary sources
nx0 = secondary_source_direction(x0);
x0 = x0(:,1:3);

% Make the size of the xs position the same as the number of secondary sources
% in order to allow x0-xs
% First we have to get the direction of a plane wave
nxs = xs / norm(xs);
xs = repmat(xs,size(x0,1),1);
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
    ls_activity = double( diag(nxs*nx0') > 0 );

elseif strcmp('ps',src) || strcmp('ls',src)
    % === Point source ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <x0-xs,n_x0> > 0
    % a = <
    %      \ 0, else
    %
    ls_activity = double( diag((x0-xs)*nx0') > 0 );

elseif strcmp('fs',src)
    % === Focused source ===
    % secondary source selection (Spors 2008)
    %
    %      / 1, if <xs-x0,n_x0> > 0
    % a = <
    %      \ 0, else
    ls_activity = double( diag((xs-x0)*nx0') > 0 );
else
    error('%s: %s is not a supported source type!',upper(mfilename),src);
end
