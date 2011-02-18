function a = wfs_amplitude(x0,y0,X,Y,xs,ys,src,conf)
%WFS_AMPLITUDE Calculates the amplitude factor for delay and add WFS simulation
%   Usage: a = wfs_amplitude(x0,y0,X,Y,xs,ys,src,conf)
%          a = wfs_amplitude(x0,y0,X,Y,xs,ys,src)
%
%   Input parameters:
%       x0, y0  - position of the loudspeaker (m)
%       X, Y    - position of the listener (m)
%       xs, ys  - position of the virtual source (m)
%       src     - source type: 'pw' - plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       conf    - optional configuration struct
%
%   Output parameters:
%       a       - amplitude factor
%
%   WFS_AMPLITUDE(x0,y0,X,Y,xs,ys,src) calculates the amplitude factor for 
%   the given loudspeaker to use for delay and add creation of a
%   loudspeaker array. 
%
%   see also: wfs_brs, echo_direction
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));
isargscalar(x0,y0,X,Y,xs,ys);
isargchar(src);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
yref = conf.yref;


%% ===== Computation =====================================================
%
%% Some thoughts about the geometry for a linear array
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
% After Spors et al. (2008) the traditional WFS formulation of the
% driving function D for a linear loudspeaker array parrallel to the 
% x-axis reads:
% D = cos(phi) / sqrt(R2),
% where phi denotes the angle between the y-axis and R2.
% Therefore we have cos(phi) = (ys-y0) / R2 and
% D = (ys-y0) .* R2.^(-3/2).

% === Amplitude factor ===
% Constant amplitude factor
g0 = sqrt(2*pi*abs(yref-y0));
if strcmp('pw',src)

    % ===== PLANE WAVE ===================================================
    % Use the position of the source as the direction vector for a plane wave
    ny0 = ys / sqrt(xs^2+ys^2);
    a = g0 * ny0;

elseif strcmp('ps',src)

    % ===== POINT SOURCE =================================================
    % Using a point source as a model we will get
    a = g0/2*pi * (y0-ys)/norm([x0 y0]-[xs ys])^2;

elseif strcmp('fs',src)

    % ===== FOCUSED SOURCE ===============================================
    % Using a amplitude corrected line sink as source model (Spors2009)
    a = g0/sqrt(2*pi) * (ys-y0)/norm([x0 y0]-[xs ys])^(3/2);

else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
end

% Add a 1/r term for the decreasing of the sound on its way from
% the loudspeaker to the listener [this comes from the 3D Green's function].
a = a * 1/norm([x0 y0]-[X Y]);
