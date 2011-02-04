function a = wfs_amplitude_linear(x0,y0,X,Y,xs,ys)
%WFS_AMPLITUDE_LINEAR Calculates the amplitude factor for delay and add WFS
%creation
%   Usage: a = wfs_amplitude_linear(x0,y0,X,Y,xs,ys)
%
%   Input parameters:
%       x0, y0  - position of the loudspeaker
%       X, Y    - position of the listener
%       xs, ys  _ position of the virtual source
%
%   Output parameters:
%       a       - amplitude factor
%
%   WFS_AMPLITUDE_LINEAR(x0,y0,X,Y,xs,ys) calculates the amplitude factor for 
%   the given loudspeaker to use for delay and add creation of a
%   linear loudspeaker array. 
%
%   see also: wfs_brs, echo_direction
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
isargscalar({x0,y0,X,Y,xs,ys},{'x0','y0','X','Y','xs','ys'});


%% ===== Computation =====================================================
% Geometry
%          [x0,y0]                [X0,Y0]
% x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^--
%             |                      |
%         R2 |  |                    |
%           |     | R                |      
%          x        |                |
%       [xs,ys]       |              |
%                       O            |
%                      [X,Y]         |
%                                    v
%                                  y-axis
%
% Distance between loudspeaker and listener: R
% Distance between loudspeaker and virtual source: R2
R = norm([X Y] - [x0 y0]);
R2 = norm([xs ys] - [x0 y0]);

% === Amplitude factor ===
% After Spors et al. (2008) the traditional WFS formulation of the
% driving function D for a linear loudspeaker array parrallel to the 
% x-axis reads:
% D = cos(phi) / sqrt(R2),
% where phi denotes the angle between the y-axis and R2.
% Therefore we have cos(phi) = (ys-y0) / R2 and
% D = (ys-y0) .* R2.^(-3/2).
% The 1/R term is for the decreasing of the sound on its way from
% the loudspeaker to the listener (R) [this comes from the 3D Green's function].
a = (ys-y0) .* R2.^(-3/2) .* 1/R;
% Amplitude for a focused source for the Driving function calculated by Hagen
%a = (ys-y0) .* R2.^(-2) .* 1/R;
