function movie_wave_field_mono_wfs_25d(X,Y,xs,L,f,src,outfile,conf)
%MOVIE_WAVE_FIELD_MONO_WFS_25D generates movie a 2.5D WFS wave field
%   Usage: movie_wave_field_mono_wfs_25d(X,Y,xs,L,f,src,outfile,conf)
%          movie_wave_field_mono_wfs_25d(X,Y,xs,L,f,src,outfile)
%
%   Input parameters:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - position of point source (m)
%       L           - array length (m)
%       f           - monochromatic frequency (Hz)
%       src         - sourcetype of the virtual source:
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       outfile     - name for the movie file
%       conf        - optional configuration struct (see SFS_config)
%
%   MOVIE_WAVE_FIELD_MONO_WFS_25D(X,Y,xs,L,f,src,outfile,conf) generates a
%   movie of simulations of a wave field of the given source positioned at xs
%   using a WFS 2.5 dimensional driving function in the temporal domain with
%   different phase.
%
%   see also: wave_field_mono_wfs_25d, plot_wavefield

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));
isargvector(X,Y);
isargposition(xs);
xs = position_vector(xs);
isargpositivescalar(L,f);
isargchar(src,outfile);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

% Plotting
useplot = conf.useplot;
% Temporary dir
tmpdir = conf.tmpdir;


%% ===== Simulation =====================================================
phase = linspace(0,2*pi,25);
% Generate a random number string for the tmp files
rn = sprintf('%04.0f',10000*rand);
% Simulate the time by different phase values
for ii = 1:length(phase)-1
    conf.phase = phase(ii);
    conf.useplot = 0;
    % Calculate wave field for the given phase
    [x,y,P] = wave_field_mono_wfs_25d(X,Y,xs,L,f,src,conf);
    x0 = secondary_source_positions(L,conf);
    ls_activity = secondary_source_selection(x0,xs,src);
    % Generate tapering window
    win = tapwin(L,ls_activity,conf);
    ls_activity = ls_activity .* win;

    % === Save temporary data ===
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    pngfile = sprintf('%s/%s_%i.png',tmpdir,rn,ii+10);
    conf.plot.mode = 'png';
    plot_wavefield(x,y,P,L,ls_activity,pngfile,conf);
end


%% ===== Create movie ====================================================
conf.useplot = useplot;
generate_movie(outfile,tmpdir,rn);

% Clean up tmp files
delete([tmpdir,'/',rn,'*.png']);


%% ===== Show movie ======================================================
if useplot
    [status,mplayer] = system('which mplayer');
    if status
        error('%s: mplayer is needed to show this movie.',upper(mfilename));
    else
        cmd = sprintf('mplayer %s -loop 0',outfile);
        system(cmd);
    end
end
