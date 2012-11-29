function movie_wave_field_imp_wfs_25d(X,Y,xs,src,L,outfile,conf)
%MOVIE_WAVE_FIELD_IMP_WFS_25D generates movie a 2.5D WFS wave field
%
%   Usage: movie_wave_field_imp_wfs_25d(X,Y,xs,src,L,outfile,[conf])
%
%   Input parameters:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - position of point source (m)
%       src         - sourcetype of the virtual source:
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       L           - array length (m)
%       outfile     - name for the movie file
%       conf        - optional configuration struct (see SFS_config)
%
%   MOVIE_WAVE_FIELD_IMP_WFS_25D(X,Y,xs,src,L,outfile,conf) generates a
%   movie of simulations of a wave field of the given source positioned at xs, ys
%   using a WFS 2.5 dimensional driving function in the temporal domain with
%   different phase.
%
%   see also: wave_field_imp_wfs_25d, plot_wavefield, generate_movie

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y);
isargposition(xs);
xs = position_vector(xs);
isargpositivescalar(L);
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
%phase = linspace(0,2*pi,25);
frame = round(linspace(-200,900,7*25));
% Generate a random number string for the tmp files
rn = sprintf('%04.0f',10000*rand);
% Disable the empty wave field warning
warning('off','SFS:check_wave_field');
% Simulate the time by different phase values
for ii = 1:length(frame)-1
    conf.frame = frame(ii);
    conf.useplot = 0;
    % Calculate wave field for the given phase
    [x,y,P] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
    x0 = secondary_source_positions(L,conf);
    x0 = secondary_source_selection(x0,xs,src);
    % Generate tapering window
    ls_activity = tapering_window(x0,conf);

    % === Save temporary data ===
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    pngfile = sprintf('%s/%s_%04.0f.png',tmpdir,rn,ii);
    conf.plot.mode = 'png';
    conf.plot.usedb = 1;
    plot_wavefield(x,y,P,L,ls_activity,pngfile,conf);
end
% Enable the empty wave field warning
warning('on','SFS:check_wave_field');

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
