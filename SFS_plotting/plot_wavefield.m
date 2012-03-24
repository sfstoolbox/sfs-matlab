function plot_wavefield(x,y,P,varargin)
%PLOT_WAVEFIELD plot the given wavefield
%
%   Usage: plot_wavefield(x,y,P,L,ls_activity,conf)
%          plot_wavefield(x,y,P,L,ls_activity)
%          plot_wavefield(x,y,P)
%
%   Input parameters:
%       x,y         - vectors for the x- and y-axis
%       P           - matrix containing the wavefield in the format P = P(y,x)
%       L           - array length. If this is given and the distance between the
%                     loudspeaker is greater than 10cm the loudspeaker are added 
%                     to the plot at their real positions.
%       ls_activity - activity of the single loudspeakers. If all loudspeakers
%                     should be active, you can simply set it to 1. Otherwise a
%                     vector with entries for every single loudspeaker is
%                     needed.
%       conf        - optional struct containing configuration variables (see
%                     SFS_config for default values)
%
%   PLOT_WAVEFIELD(x,y,P,L,ls_activity) plots the wavefield P in dependence of 
%   the x and y axes. Therefore the wavefield is normalized to 1 at its center
%   position P(end/2,end/2). For a given array length L also the loudspeaker are
%   added to the plot at their real positions. But only if distance between them
%   is larger than 10cm.
%
%   see also: 

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));
isargvector(x,y);
isargmatrix(P);
% Setting default values and check varargin parameters
ls_activity = 1;
% FIXME: remove outfile here, because we now use conf.plot.file
outfile = 'wavefield.plot';
for ii = 1:length(varargin)
    if isstruct(varargin{ii})
        conf = varargin{ii};
    elseif ischar(varargin{ii})
        outfile = varargin{ii};
    elseif isscalar(varargin{ii}) && varargin{ii}>0 && ~exist('L','var')
        L = varargin{ii};
        if ii+1<=length(varargin) && isvector(varargin{ii})
            ls_activity = varargin{ii+1};
        end
    end
end
if ~exist('conf','var')
    conf = SFS_config;
end


%% ===== Configuration ==================================================

% Check if we have a loudspeaker array at all
if ~exist('L','var')
    conf.plot.loudspeakers = 0;
end

% Tmp dir  
tmpdir = conf.tmpdir;
% Center position of array
X0 = position_vector(conf.X0);
% Distance between loudspeakers
dx0 = conf.dx0;
% Used array geometry
array = conf.array;
% Plotting
useplot = conf.useplot;
p.usegnuplot = conf.plot.usegnuplot;
p.cmd = conf.plot.cmd;
p.usedb = conf.plot.usedb;
p.mode = conf.plot.mode;
p.size = conf.plot.size;
p.caxis = conf.plot.caxis;
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

% Check if the array length is given, if the loudspeaker should be plotted
if ~exist('L','var') && p.loudspeakers
    error(['%s: the array length L has to be given in order to draw ',...
           'the loudspeakers!'],upper(mfilename));
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
    GraphDefaults(p.mode);

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

    % Change colormap to Gray
    % NOTE: Earlier versions of Octave didn't know the Gray colormap.
    % So we have to set the values by hand.
    % The values are generated with map = colormap('Gray') in Matlab.
    %colormap(gray);
    map = [0  0         0
    0.0159    0.0159    0.0159
    0.0317    0.0317    0.0317
    0.0476    0.0476    0.0476
    0.0635    0.0635    0.0635
    0.0794    0.0794    0.0794
    0.0952    0.0952    0.0952
    0.1111    0.1111    0.1111
    0.1270    0.1270    0.1270
    0.1429    0.1429    0.1429
    0.1587    0.1587    0.1587
    0.1746    0.1746    0.1746
    0.1905    0.1905    0.1905
    0.2063    0.2063    0.2063
    0.2222    0.2222    0.2222
    0.2381    0.2381    0.2381
    0.2540    0.2540    0.2540
    0.2698    0.2698    0.2698
    0.2857    0.2857    0.2857
    0.3016    0.3016    0.3016
    0.3175    0.3175    0.3175
    0.3333    0.3333    0.3333
    0.3492    0.3492    0.3492
    0.3651    0.3651    0.3651
    0.3810    0.3810    0.3810
    0.3968    0.3968    0.3968
    0.4127    0.4127    0.4127
    0.4286    0.4286    0.4286
    0.4444    0.4444    0.4444
    0.4603    0.4603    0.4603
    0.4762    0.4762    0.4762
    0.4921    0.4921    0.4921
    0.5079    0.5079    0.5079
    0.5238    0.5238    0.5238
    0.5397    0.5397    0.5397
    0.5556    0.5556    0.5556
    0.5714    0.5714    0.5714
    0.5873    0.5873    0.5873
    0.6032    0.6032    0.6032
    0.6190    0.6190    0.6190
    0.6349    0.6349    0.6349
    0.6508    0.6508    0.6508
    0.6667    0.6667    0.6667
    0.6825    0.6825    0.6825
    0.6984    0.6984    0.6984
    0.7143    0.7143    0.7143
    0.7302    0.7302    0.7302
    0.7460    0.7460    0.7460
    0.7619    0.7619    0.7619
    0.7778    0.7778    0.7778
    0.7937    0.7937    0.7937
    0.8095    0.8095    0.8095
    0.8254    0.8254    0.8254
    0.8413    0.8413    0.8413
    0.8571    0.8571    0.8571
    0.8730    0.8730    0.8730
    0.8889    0.8889    0.8889
    0.9048    0.9048    0.9048
    0.9206    0.9206    0.9206
    0.9365    0.9365    0.9365
    0.9524    0.9524    0.9524
    0.9683    0.9683    0.9683
    0.9841    0.9841    0.9841
    1.0000    1.0000    1.0000];
    % Invert white and black color
    map = map(end:-1:1,:,:);
    colormap(map);

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
            x0 = secondary_source_positions(L,conf);
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
        x0 = secondary_source_positions(L,conf);

        % FIXME: I think the following code is not neccessary
        % check if we have loudspeakers outside of the desired plotting region and
        % remove them in order to don't mess up the graph
        %idx = (( x0<x(1) || x0>x(end) || y0<y(1) || y0>y(end) ));
        %x0(idx) = NaN;
        %y0(idx) = NaN;
        %phi(idx) = NaN;

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
            phi = cart2pol(x0(:,4)-x0(:,1),x0(:,5)-x0(:,2));
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

    elseif strcmp('monitor',p.mode) | strcmp('screen',p.mode)
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
