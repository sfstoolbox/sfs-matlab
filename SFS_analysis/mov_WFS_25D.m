function [x,y] = mov_WFS_25D(X,Y,xs,ys,L,f,src,outfile,conf)
%MOV_WFS_25D generates a wave field movie for the 2.5D case using WFS
%   Usage: mov_WFS_25D(X,Y,xs,ys,L,f,src,outfile,conf)
%          mov_WFS_25D(X,Y,xs,ys,L,f,src,outfile)
%
%   Input parameters:
%       X           - length of the X axis (m) [xaxis: -X/2:X/2]
%       Y           - length of the Y axis (m) [yaxis: -0.1:Y]
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
%       L           - array length (m)
%       f           - Monochromatic frequency (Hz)
%       src         - sourcetype of the virtual source:
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       outfile     - name for the movie file
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%
%   MOV_WFS_25D(X,Y,xs,ys,L,f,src,outfile,conf) generates a movie of
%   simulations of a wave field of the given source positioned at xs, ys
%   using a WFS 2.5 dimensional driving function in the temporal domain with 
%   different phase.
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%
%   see also: wf_WFS_25D, plot_wavefield

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 8;
nargmax = 9;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(X) || ~isvector(X)
    error('%s: X has to be a vector!',upper(mfilename));
end
if ~isnumeric(Y) || ~isvector(Y)
    error('%s: Y has to be a vector!',upper(mfilename));
end
if ~isnumeric(xs) || ~isscalar(xs)
    error('%s: xs has to be a scalar!',upper(mfilename));
end
if ~isnumeric(ys) || ~isscalar(ys)
    error('%s: ys has to be a scalar!',upper(mfilename));
end
if ~isnumeric(L) || ~isscalar(L) || L<=0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end
if ~isnumeric(f) || ~isscalar(f) || f<=0
    error('%s: f has to be a positive scalar!',upper(mfilename));
end
if ~ischar(src)
    error('%s: src has to be a string!',upper(mfilename));
end
if ~ischar(outfile)
    error('%s: outfile has to be a string!',upper(mfilename));
end
if nargin<nargmax
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
if(useconfig)
    conf = SFS_config;
end

% Array position (m)
X0 = conf.X0;                    
Y0 = conf.Y0;

% SFS Path
sfspath = conf.sfspath;

% Plotting
useplot = conf.useplot;
usegnuplot = conf.usegnuplot;

% Temporary dir
tmpdir = conf.tmpdir;


%% ===== Simulation =====================================================
phase = linspace(0,2*pi,25);
% Simulate the time by different phase values
for i = 1:length(phase)-1
    conf.phase = phase(i);
    conf.useplot = 0;
    % Calculate wave field for the given phase
    [x,y,P] = wf_WFS_25D(X,Y,xs,ys,L,f,src,conf);
    % Replace the wave field with zeros for y<0
    % Find y<0
    idx = (( y<0 ));
    % Set wave field to zero in this area
    P(:,idx) = 0;

    % === Save data ===
    mkdir(tmpdir);
    datafile = sprintf('%s/%i.dat',tmpdir,i+10);
    if(usegnuplot)
        gp_save_matrix(datafile,x,y,real(P))
    else
        matlab_version_missing(mfilename);
    end
end


%% ===== Create movie ====================================================
conf.useplot = useplot;
movie_wavefield(x,y,L,outfile,conf);

% Clean up tmp files
delete([tmpdir,'/*']);
rmdir(tmpdir);
