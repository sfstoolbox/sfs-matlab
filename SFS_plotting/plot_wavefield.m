function plot_wavefield(x,y,z,P,x0,ls_activity,conf)
%PLOT_WAVEFIELD plot the given wavefield
%
%   Usage: plot_wavefield(x,y,z,P,[x0,[ls_activity]],[conf])
%
%   Input parameters:
%       x,y,z       - vectors for the x-, y- and z-axis
%       P           - matrix containing the wavefield in the format P = P(y,x)
%       x0          - matrix containing the secondary source positions to plot.
%                     Default: plot no secondary sources
%       ls_activity - activity of the single secondary sources. Note: this option
%                     is only valid, if you give also x0 as an input parameter.
%                     The default behavior is to plot all secondary sources as 
%                     active.
%       conf        - optional configuration struct (see SFS_config)
%
%   PLOT_WAVEFIELD(x,y,z,P,L,ls_activity,conf) plots the wavefield P in dependence
%   of the x and y axes. Therefore the wavefield is normalized to 1 at its
%   center position P(end/2,end/2). For a given set x0 of secondary sources the
%   loudspeakers are added to the plot at their real positions. But only if
%   distance between them is larger than 10cm. The ls_activity option specifies
%   the color shade of the speakers, going from 0 (white) to 1(dark gray). The
%   default behavior is to set all speakers to 1.
%
%   see also: wave_field_mono_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
nargmin = 4;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(x,y,z);
isargmatrix(P);
if nargin==nargmax-1
    if isstruct(ls_activity)
        conf = ls_activity;
        ls_activity = ones(1,size(x0,1));
    else
        conf = SFS_config;
    end
elseif nargin==nargmax-2
    if isstruct(x0)
        conf = x0;
        x0 = [];
    else
        conf = SFS_config;
        ls_activity = ones(1,size(x0,1));
    end
elseif nargin==nargmax-3
    conf = SFS_config;
end
if ~exist('x0','var') || length(x0)==0
    conf.plot.loudspeakers = 0;
else
%      isargsecondarysource(x0);
end
isargstruct(conf);


%% ===== Configuration ==================================================
dx0 = conf.dx0;
% Tmp dir
tmpdir = conf.tmpdir;
% Center position of array
X0 = conf.X0;
% Plotting
p.usegnuplot = conf.plot.usegnuplot;
p.cmd = conf.plot.cmd;
p.usedb = conf.plot.usedb;
p.mode = conf.plot.mode;
p.size = conf.plot.size;
p.size_unit = conf.plot.size_unit;
p.caxis = conf.plot.caxis;
p.colormap = conf.plot.colormap;
p.loudspeakers = conf.plot.loudspeakers;
p.realloudspeakers = conf.plot.realloudspeakers;
p.lssize = conf.plot.lssize;
p.usefile = conf.plot.usefile;
p.file = conf.plot.file;


%% ===== Calculation =====================================================
% Handle the given axis and check which should be plotted
[dimensions,x1,x2] = xyz_axes_selection(x,y,z);
if ~dimensions(1)
    % FIXME: in order to work with gnuplot the label should be prtinted
    % with the extra function, which can handle if the output should be
    % LaTeX or something else
    %str_xlabel = print_label('y','m',conf);
    p.xlabel = 'y / m';
    p.ylabel = 'z / m';
elseif ~dimensions(2)
    p.xlabel = 'x / m';
    p.ylabel = 'z / m';
elseif ~dimensions(3)
    p.xlabel = 'x / m';
    p.ylabel = 'y / m';
else
    % FIXME: in this case every three axis should be plotted and we should
    % switch to use splot or some other alternativ to plot it in 3D.
    to_be_implemented(mfilename);
end

% Check the size of x,y and P
if size(P,1)~=length(x2) || size(P,2)~=length(x1)
    error('%s: the size of P has to be x2 x x1.',upper(mfilename));
end

if(p.usedb)
    % Check if we have any activity in the wave field
    if max(abs(P(:)))~=0
        % For the dB case scale the signal maximum to 0 dB
        %P = P./max(abs(P(:)));
    else
        % If we have only zeros in the wave field set the field to eps to avoid
        % problems with log(0).
        P(:) = eps;
    end
end

% Check if we should plot loudspeakers symbols and fix the size of ls_activity
if p.loudspeakers && dx0<=0.01
    warning(['%s: the given loudspeaker distance is to small. ',...
            'Disabling plotting of the loudspeakers'],upper(mfilename));
    p.loudspeakers = 0;
else
    % fixing the length of ls_activity
    if length(ls_activity)==1
        ls_activity = repmat(ls_activity,[1 size(x0,1)]);
    end
end

% set the color bar axis to default values if not given otherwise
if p.caxis else
    if p.usedb
        p.caxis = [-45,0];
    else
        p.caxis = [-1,1];
    end
end


%% ===== Plotting ========================================================

if ~(p.usegnuplot)
    % ===== Plot the wave field with Matlab/Octave =======================
    %
    % Create a new figure
    figure;
    % set size
    figsize(p.size(1),p.size(2),p.size_unit);

    % Plotting
    if(p.usedb)
        % Plot the amplitude of the wave field in dB
        imagesc(x1,x2,20*log10(abs(P)),p.caxis);
    else
        % Plot the wave field
        imagesc(x1,x2,real(P),p.caxis);
    end

    % Add color bar
    set_colorbar(conf);

    % Set the y direction in normal mode (imagesc uses the reverse mode by
    % default)
    turn_imagesc;

    % Set the axis to use the same amount of space for the same length (m)
    axis image;
    % Labels etc. for the plot
    xlabel(p.xlabel);
    ylabel(p.ylabel);

    % Add loudspeaker to the plot
    if p.loudspeakers && dimensions(1) && dimensions(2)
        hold on;
        draw_loudspeakers(x0,ls_activity,conf);
        hold off;
    end

    % Save as file
    if ~isempty(p.file) && strcmp('png',p.file(end-2:end))
        print_png(p.file,conf);
    elseif ~isempty(p.file) && strcmp('eps',p.file(end-2:end))
        print_eps(p.file,conf);
    end

else


%% ===== Plot the wave field using Gnuplot ===============================

    % ----- Store all the files needed for plotting ----------------------
    % tmp dir for storing temporary files
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    % Create output file name
    if p.usefile
        p.datafile = sprintf('%s.dat',p.file(end-2:end));
        p.lsfile = sprintf('%s_ls.txt',p.file(end-2:end));
        p.gnuplotfile = sprintf('%s.gnu',p.file(end-2:end));
    else
        % Generate a random number string for the tmp files
        rn = sprintf('%04.0f',10000*rand);
        p.datafile = sprintf('%s/wavefield%s.dat',tmpdir,rn);
        p.lsfile = sprintf('%s/loudspeakers%s.txt',tmpdir,rn);
        p.gnuplotfile = sprintf('%s/gnuplot%s.gnu',tmpdir,rn);
    end
    % Storing loudspeaker positions and activity
    if(p.loudspeakers)
        [phi,~] = cart2pol(x0(:,4),x0(:,5));
        [x0,y0,phi,ls_activity] = column_vector(x0(:,1),x0(:,2),phi,ls_activity);
        gp_save(p.lsfile,x0,[y0 phi ls_activity]);
    end

    % Check if we should handle the wave field in dB
    if p.usedb
        % Save the data for plotting with Gnuplot
        gp_save_matrix(p.datafile,x1,x2,db(abs(P)));
        p.cbtics = 5;
        p.dim = 'p';
        p.unit = 'dB';
    else
        % Save the data for plotting with Gnuplot
        gp_save_matrix(p.datafile,x1,x2,real(P));
        p.cbtics = 1;
        p.dim = 'P';
        p.unit = '';
    end

    % stor the x- and y-range
    p.xmin = x1(1);
    p.xmax = x1(end);
    p.ymin = x2(1);
    p.ymax = x2(end);
    % plot the files
    if ~isempty(p.file) && strcmp('png',p.file(end-2:end))
        % png file
        fprintf(1,'Plotted to file %s\n',p.file);
        cmd = gp_print_png(p);
    elseif ~isempty(p.file) && strcmp('eps',p.file(end-2:end))
        % eps file
        fprintf(1,'Plotted to file %s\n',p.file);
        cmd = gp_print_eps(p);
    elseif ~isempty(p.file) && strcmp('tex',p.file(end-2:end))
        % epslatex file
        fprintf(1,'Plotted to file %s\n',p.file);
        cmd = gp_print_epslatex(p);
    else
        % on ths screen
        cmd = gp_print_screen(p);
    end


    if p.usefile
        fid = fopen(p.gnuplotfile,'w');
        fprintf(fid,'%s',cmd);
        fclose(fid);
        run_cmd = sprintf('gnuplot %s\n',p.gnuplotfile);
    else
        run_cmd = sprintf('gnuplot<<EOC\n%s\nEOC\n',cmd);
    end
    % Start Gnuplot for plotting the data
    system(run_cmd);

    if ~p.usefile
        % Remove tmp files
        if exist(p.datafile,'file')
            delete(p.datafile);
        end
        if exist(p.lsfile,'file')
            delete(p.lsfile);
        end
        if exist(p.gnuplotfile,'file')
            delete(p.gnuplotfile);
        end
    end

end
