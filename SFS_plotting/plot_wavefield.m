function plot_wavefield(x,y,P,x0,ls_activity,conf)
%PLOT_WAVEFIELD plot the given wavefield
%
%   Usage: plot_wavefield(x,y,P,[x0,[ls_activity]],[conf])
%
%   Input parameters:
%       x,y         - vectors for the x- and y-axis
%       P           - matrix containing the wavefield in the format P = P(y,x)
%       x0          - matrix containing the secondary source positions to plot.
%                     Default: plot no secondary sources
%       ls_activity - activity of the single secondary sources. Note: this option
%                     is only valid, if you give also x0 as an input parameter.
%                     The default behavior is to plot all secondary sources as 
%                     active.
%       conf        - optional configuration struct (see SFS_config)
%
%   PLOT_WAVEFIELD(x,y,P,x0,ls_activity,conf) plots the wavefield P in dependence
%   of the x and y axes. Therefore the wavefield is normalized to 1 at its
%   center position P(end/2,end/2). For a given set x0 of secondary sources the
%   loudspeakers are added to the plot at their real positions. But only if
%   distance between them is larger than 10cm. The ls_activity option specifies
%   the color shade of the speakers, going from 0 (white) to 1(dark gray). The
%   default behavior is to set all speakers to 1.
%
%   see also: wave_field_mono_wfs_25d

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
nargmin = 3;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(x,y);
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
    end
elseif nargin==nargmax-3
    conf = SFS_config;
end
if ~exist('x0','var') || length(x0)==0
    conf.plot.loudspeakers = 0;
else
    isargsecondarysource(x0);
end
isargstruct(conf);


%% ===== Configuration ==================================================
dx0 = conf.dx0;
% Tmp dir
tmpdir = conf.tmpdir;
% Center position of array
X0 = position_vector(conf.X0);
% Plotting
useplot = conf.useplot;
p.usegnuplot = conf.plot.usegnuplot;
p.cmd = conf.plot.cmd;
p.usedb = conf.plot.usedb;
p.mode = conf.plot.mode;
p.size = conf.plot.size;
p.caxis = conf.plot.caxis;
p.colormap = conf.plot.colormap;
p.loudspeakers = conf.plot.loudspeakers;
p.realloudspeakers = conf.plot.realloudspeakers;
p.lssize = conf.plot.lssize;
p.usefile = conf.plot.usefile;
p.file = conf.plot.file;


%% ===== Calculation =====================================================
% Check the size of x,y and P
if size(P,1)~=length(y) || size(P,2)~=length(x)
    error('%s: the size of P has to be y x x.',upper(mfilename));
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


%% ===== Plotting ========================================================

if ~(p.usegnuplot)
    % ===== Plot the wave field with Matlab/Octave =======================
    %
    % Create a new figure
    figure;

    if(p.usedb)
        % Plot the amplitude of the wave field in dB
        imagesc(x,y,20*log10(abs(P)),[-45 0]);
        % Set the limits of the colormap and add a colorbar
        if length(p.caxis)==2
            caxis(p.caxis);
        end
        h = colorbar;
        ylabel(h,'Amplitude (dB)');
        % Get the font size and name of the figure and adjust the colorbar
        fsize = get(gca,'FontSize');
        fname = get(gca,'FontName');
        set(h,'FontSize',fsize);
        set(h,'FontName',fname);
        temp = get(h,'Ylabel');
        set(temp,'FontSize',fsize);
        set(temp,'FontName',fname);
    else
        % Plot the wave field
        imagesc(x,y,real(P),[-1 1]);
        % Set the limits of the colormap and add a colorbar
        if length(p.caxis)==2
            caxis(p.caxis);
        end
        colorbar;
    end

    % Set the y direction in normal mode (imagesc uses the reverse mode by
    % default)
    turn_imagesc;

    % Change colormap (default: gray)
    if p.colormap
        colormap(p.colormap);
    else
        colormap('gray');
    end

    % Set the axis to use the same amount of space for the same length (m)
    axis image;
    % Labels etc. for the plot
    xlabel('x (m)');
    ylabel('y (m)');

    % Add loudspeaker to the plot
    if(p.loudspeakers)
        if dx0<=0.01
            warning(['%s: the given loudspeaker distance is to small. ',...
                     'Disabling plotting of the loudspeakers'],upper(mfilename));
        else
            hold on;
            draw_loudspeakers(x0,ls_activity,conf);
            hold off;
        end
    end

    if strcmp('png',p.mode)
        print(outfile,'-dpng','-S640,480');
        close;
    elseif strcmp('paper',p.mode)
        print(outfile,'-deps','-S320,240');
        close;
    end

else


%% ===== Plot the wave field using Gnuplot ===============================

    % tmp dir for storing temporary files
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end

    % Create output file name
    if p.usefile
        datafile = sprintf('%s.dat',p.file);
        lsfile = sprintf('%s_ls.txt',p.file);
    else
        % Generate a random number string for the tmp files
        rn = sprintf('%04.0f',10000*rand);
        datafile = sprintf('%s/wavefield%s.dat',tmpdir,rn);
        lsfile = sprintf('%s/loudspeakers%s.txt',tmpdir,rn);
    end

    % Check if we should plot the loudspeakers.
    if(p.loudspeakers)
        % Loudspeaker positions and directions
        if  dx0<= 0.01
            warning(['%s: the given loudspeaker distance is to small. ',...
                    'Disabling plotting of the loudspeakers'],upper(mfilename));
            p.loudspeakers = 0;
        else
            % fixing the length of ls_activity
            if length(ls_activity)==1
                ls_activity = repmat(ls_activity,size(x0));
            end
            % Storing loudspeaker positions and activity
            [phi,r_tmp] = cart2pol(x0(:,4),x0(:,5));
            [x0,y0,phi,ls_activity] = column_vector(x0(:,1),x0(:,2),phi,ls_activity);
            gp_save(lsfile,x0,[y0 phi ls_activity]);
        end
    end

    % Check if we should handle the wave field in dB
    if p.usedb
        % Save the data for plotting with Gnuplot
        gp_save_matrix(datafile,x,y,db(abs(P)));
        if p.caxis else
            p.caxis = [-45,0];
        end
        cbtics = 5;
        pdim = 'p';
        punit = 'dB';
    else
        % Save the data for plotting with Gnuplot
        gp_save_matrix(datafile,x,y,real(P));
        if p.caxis else
            p.caxis = [-1,1];
        end
        cbtics = 1;
        pdim = 'P';
        punit = '';
    end

    %% === set common Gnuplot commands
    cmd = sprintf([...
        '#!/usr/bin/gnuplot\n', ...
        '# generated by plot_wavefield.m\n', ...
        'unset key\n', ...
        'set size ratio -1\n\n', ...
        '# border\n', ...
        'set style line 101 lc rgb ''#808080'' lt 1 lw 1\n', ...
        'set border front ls 101\n\n', ...
        'set colorbox\n', ...
        'set palette gray negative\n', ...
        'set xrange [%f:%f]\n', ...
        'set yrange [%f:%f]\n', ...
        'set cbrange [%f:%f]\n', ...
        'set tics scale 0.75\n', ...
        'set cbtics scale 0\n', ...
        'set xtics 1\n', ...
        'set ytics 1\n', ...
        'set cbtics %f\n', ...
        'set xlabel ''%s''\n', ...
        'set ylabel ''%s''\n', ...
        'set label ''%s'' at screen 0.84,0.14\n'], ...
        x(1),x(end), ...
        y(1),y(end), ...
        p.caxis(1),p.caxis(2), ...
        cbtics, ...
        print_label('x','m',conf), ...
        print_label('y','m',conf), ...
        print_label(pdim,punit,conf));


    if strcmp('paper',p.mode)
        %% === Paper ===
        % Therefore we use the epslatex terminal of Gnuplot, see:
        % http://www.gnuplotting.org/introduction/output-terminals/#epslatex
        cmd = sprintf([...
            '%s\n', ...
            'set t epslatex size %fcm,%fcm color colortext\n', ...
            'set output ''%s.tex'';\n', ...
            'set style line 1 lc rgb ''#000000'' pt 2 ps 2 lw 2\n', ...
            'set format ''$%%g$''\n\n', ...
            '%s\n', ...
            '%s\n'], ...
            cmd, ...
            p.size(1),p.size(2), ...
            p.file, ...
            p.cmd);

    elseif strcmp('talk',p.mode)
        %% === Talk ===
        % Therefore we use the epslatex terminal of Gnuplot, see:
        % http://www.gnuplotting.org/introduction/output-terminals/#epslatex
        cmd = sprintf([...
            '%s\n', ...
            'set t epslatex size %fcm,%fcm color colortext\n', ...
            'set output ''%s.tex''\n', ...
            'set style line 1 lc rgb ''#000000'' pt 2 ps 2 lw 3\n', ...
            'set format ''$%%g$''\n\n', ...
            '%s\n', ...
            '%s\n'], ...
            cmd, ...
            p.size(1),p.size(2), ...
            p.file, ...
            p.cmd);

    elseif strcmp('monitor',p.mode) || strcmp('screen',p.mode)
        % === Monitor ===
        % Therefore we use the wxt Gnuplot terminal
        cmd = sprintf([...
            '%s\n', ...
            'set t wxt size 700,524 enhanced font ''Verdana,14'' persist\n', ...
            'set style line 1 lc rgb ''#000000'' pt 2 ps 2 lw 2\n\n', ...
            '%s\n'], ...
            cmd, ...
            p.cmd);

    elseif strcmp('png',p.mode)
        % === png ====
        % Therefore we use the pngcairo Gnuplot terminal, see:
        % http://www.gnuplotting.org/introduction/output-terminals/#pngsvg
        cmd = sprintf([...
            '%s\n', ...
            'set t pngcairo size %fcm,%fcm enhanced font ''Verdana,12'';\n', ...
            'set output ''%s.png'';\n', ...
            'set style line 1 lc rgb ''#000000'' pt 2 ps 2 lw 2;\n\n', ...
            '%s\n'], ...
            cmd, ...
            p.size(1),p.size(2), ...
            p.file, ...
            p.cmd);
    else
        error('%s: %s is not a valid plotting mode!',upper(mfilename),p.mode);
    end

    % Adding loudspeaker drawing and plotting of the wave field
    if p.loudspeakers && p.realloudspeakers
        % Plotting real loudspeaker symbols
        cmd = sprintf(['%s', ...
            'call ''gp_draw_loudspeakers.gnu'' ''%s'' ''%f''\n', ...
            'plot ''%s'' binary matrix with image'], ...
            cmd,lsfile,p.lssize,datafile);
    elseif p.loudspeakers
        % Plotting only points at the loudspeaker positions
        cmd = sprintf(['%s', ...
            'plot ''%s'' binary matrix with image,\n', ...
            '     ''%s'' u 1:2 w points ls 1'], ...
            cmd,datafile,lsfile);
    else
        % plotting no loudspeakers at all
        cmd = sprintf(['%s', ...
            'plot ''%s'' binary matrix with image'], ...
            cmd,datafile);
    end

    if p.usefile
        gnufile = sprintf('%s.gnu',p.file);
        fid = fopen(gnufile,'w');
        fprintf(fid,'%s',cmd);
        fclose(fid);
    else
        cmd = sprintf('gnuplot<<EOC\n%s\nEOC\n',cmd);
        % Start Gnuplot for plotting the data
        system(cmd);
        % Remove tmp files
        if exist(datafile,'file')
            delete(datafile);
        end
        if exist(lsfile,'file')
            delete(lsfile);
        end
    end

end
