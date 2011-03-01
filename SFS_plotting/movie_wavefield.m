function movie_wavefield(x,y,L,outfile,conf)
%PLOT_WAVEFIELD plot the given wavefield
%   Usage: plot_wavefield(x,y,L,outfile,conf)
%          plot_wavefield(x,y,L,outfile)
%
%   Input parameters:
%       x,y     - vectors for the x- and y-axis
%       P       - matrix containing the wavefield
%       outfile - base name of the input files
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   MOVIE_WAVEFIELD(x,y,L,outfile) plots a wavefield movie for the field
%   P stored in the files basename[11-35] in the tmpdir directory in
%   dependence of the x and y axes. Therefore the wavefield is normalized 
%   to 1 at its center position P(end/2,end/2).
%
%   see also: plot_wavefield, WFS_25D_wavefield, SDM_25D_wavefield
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(x) || ~isvector(x)
    error('%s: x has to be a vector!',upper(mfilename));
end
if ~isnumeric(y) || ~isvector(y)
    error('%s: y has to be a vector!',upper(mfilename));
end
if ~isnumeric(L) || ~isscalar(L) || L<0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end
if ~ischar(outfile)
    error('%s: outfile has to be a string!',upper(mfilename))
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

% SFS Path
sfspath = conf.sfspath;

% Center position of array
X0 = conf.X0;
Y0 = conf.Y0;
% Distance between loudspeakers
LSdist = conf.LSdist;

% Use gnuplot?
usegnuplot = conf.usegnuplot;
usedb = conf.usedb;
useeps = conf.useeps;

useplot = conf.useplot;

% Temporary dir
tmpdir = conf.tmpdir;


%% ===== Calculation ====================================================

% Loudspeaker positions (LSdir describes the directions of the LS) for a 
% linear WFS array
nLS = fix(L/LSdist)+1;
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*LSdist,nLS);
[x0,y0,phi] = secondary_source_positions(L,conf);


%% ===== Plotting ========================================================

if(usegnuplot)

    % Store the wave field and the loudspeaker positions

    % Generate the Gnuplot command line
    cmd = sprintf(['gnuplot<<EOC\n', ...
        'set loadpath "%s/SFS_plotting"\n', ...
        'call "gp_wavefield_movie.gnu" "%f" "%f" "%f" "%f" ', ...
        '"%f" "%f" "%f" "%s"\nEOC\n'], ...
        sfspath,x(1),x(end),y(1),y(end),nLS,LSdist,LSpos(1),tmpdir);

    % Start Gnuplot for plotting the data
    system(cmd);

    % Generate a movie from the png files with MEncoder
    cmd = sprintf(['mencoder mf://%s/*.png -mf w=700:h=524:fps=25:type=png', ...
        ' -ovc copy -oac copy -o %s.avi'],tmpdir,outfile);
    system(cmd);

else
    matlab_version_missing(mfilename);
end

% Show movie
if(useplot)
    cmd = sprintf('mplayer %s.avi -loop 0',outfile);
    system(cmd);
end
