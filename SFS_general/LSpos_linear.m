function [LSpos,LSdir] = LSpos_linear(X0,Y0,size,nLS)
%LSPOS_LINEAR Generates the loudspeaker positions and directions for a
%linear WFS array
%   Usage: [LSpos,LSdir] = LSpos_linear(X0,Y0,size,nLS)
%
%   Options:
%       X0, Y0      - Position of the array center (m)
%       size        - Size of the array (m)
%       nLS         - Number of loudspeaker to use
%
%   Output:
%       LSpos       - Position of the loudspeaker (m)
%                     (LSpos(1,:) = LSpos_x, LSpos(2,:) = LSpos_y)
%       LSdir       - Directions of the loudspeakers (rad)
%                     (pi/2 for every single loudspeaker)
%
%   LSPOS_LINEAR(X0,Y0,size,nLS) generates the loudspeaker
%   positions LSpos and their directions LSdir for a linear WFS array
%
%   Geometry:     
%    x-axis                   
%       <--------------------------------------------------------
%                                   |
%          LS (Loudspeaker)         | 
%          |           [X0 Y0]      |
%          ^--^--^--^--^--^--^--^--^|-^--^                    
%                   (Array Center)  |
%                                   |
%                                   |
%                                   v y-axis
%
% see also: wfs_brs, LSpos_circle
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

if nargchk(4,4,nargin)
    error(['Wrong number of args.',...
           'Usage: [LSpos,LSdir] = LSpos_linear(X0,Y0,size,nLS)']);
end

if ~isnumeric(X0) || ~isscalar(X0)
    error('%s: X0 has to be a scalar!');
end

if ~isnumeric(Y0) || ~isscalar(Y0)
    error('%s: Y0 has to be a scalar!');
end

if ~isnumeric(size) || ~isscalar(size) || size<0
    error('%s: size has to be a positive scalar!');
end

if ~isnumeric(nLS) || ~isscalar(nLS) || nLS<0
    error('%s: nLS has to be a positive scalar');
end


%% ===== Calculation ====================================================
% Positions of the loudspeakers
LSpos(1,:) = X0 + linspace(-size/2,size/2,nLS);
% === Add jitter to the loudspeaker positions ===
%LSpos(1,1) = X0-size/2;
%LSpos(1,nLS) = X0+size/2;
%for ii=2:nLS-1
%    jitter = size/(4*nLS) * randn;
%    LSpos(1,ii) = X0-size/2+(ii-1)*size/nLS + jitter;
%end
LSpos(2,:) = Y0 * ones(1,nLS);
% Direction (orientation) of the loudspeaker
LSdir = pi*ones(1,nLS);


