function [] = draw_head(X,Y,phi)
%DRAW_HEAD draws a head symbol at the given position
%
%   Usage: draw_head(X,Y,phi)
%
%   Input options:
%       X,Y         - positions of the head (m)
%       phi         - directions of the head (rad)
%
%   DRAW_HEAD(X,Y,phi) draws a head symbols at the given position. The head
%   symbol is pointing in the direction given by phi.
%
% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameter =====================================
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargscalar(X,Y,phi);


%% ===== Plotting ========================================================
to_be_implemented;
