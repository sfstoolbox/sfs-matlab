function nx0 = secondary_source_direction(x0)
%SECONDARY_SOURCE_DIRECTION returns the direction of secondary sources
%
%   Usage: secondary_source_direction(x0)
%
%   Input parameters:
%       x0      - secondary source positions and directions (m)
%
%   SECONDARY_SOURCE_DIRECTION(x0) calculates the direction the secondary
%   sources are pointing to by substracting x0(:,4:6) and x0(:,1:3) and
%   normalize the output.
%
%   see also: secondary_source_positions, secondary_source_selection,
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin)),
isargsecondarysource(x0);


%% ===== Calculation ====================================================
% Direction of secondary sources
% FIXME: find another solution than for loop
for ii = 1:size(x0,1)
    nx0(ii,:) = (x0(ii,4:6)-x0(ii,1:3)) / norm(x0(ii,4:6)-x0(ii,1:3));
end
