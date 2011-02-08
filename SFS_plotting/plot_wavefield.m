function plot_wavefield(x,y,P,L,ls_activity,conf)
%PLOT_WAVEFIELD plot the given wavefield
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
%   see also: wf_WFS_25D, wf_SDM_25D
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
isargvector(x,y);
isargmatrix(P);
if exist('L','var')
    if isstruct(L)
        conf = L;
        clear L
    else
        isargpositivescalar(L);
    end
end
if exist('ls_activity','var')
    if isstruct(ls_activity)
        conf = ls_activity;
        ls_activity = 1;
    else
        isargvector(ls_activity);
    end
else
    ls_activity = 1;
end
if ~exist('conf','var')
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

% Check if we have a loudspeaker array at all
if ~exist('L','var')
    conf.plot.loudspeakers = 0;
end
% SFS Path
sfspath = conf.sfspath;
% Center position of array
X0 = conf.X0;
Y0 = conf.Y0;
% Distance between loudspeakers
LSdist = conf.LSdist;
% Used array geometry
array = conf.array;
% Plotting
useplot = conf.useplot;
p.usegnuplot = conf.plot.usegnuplot;
p.usedb = conf.plot.usedb;
p.mode = conf.plot.mode;
p.outfile = conf.plot.outfile;
p.caxis = conf.plot.caxis;
p.loudspeakers = conf.plot.loudspeakers;


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

% FIXME: this has also to be checked for the other array geometries!
if strcmp('linear',array)
    % Replace the wave field with zeros for y<0
    % Find y<0
    idx = (( y<0 ));
    % Set wave field to zero in this area
    P(idx,:) = 0;
end


%% ===== Plotting ========================================================

if ~(p.usegnuplot)
    % === Plot the wave field with Matlab/Octave ===
    %
    % Create a new figure
    figure;
    GraphDefaults(p.mode);

    if(p.usedb)
        % Plot the amplitude of the wave field in dB
        imagesc(x,y,20*log10(abs(P)),[-35 5]);
        % Set the limits of the colormap and add a colorbar
        caxis(p.caxis);
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
        caxis(p.caxis);
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
        if LSdist<=0.01
            warning(['%s: the given loudspeaker distance is to small. ',...
                     'Disabling plotting of the loudspeakers'],upper(mfilename));
        else
            [x0,y0,phi] = secondary_source_positions(L,conf);
            hold on;
            draw_loudspeakers(x0,y0,phi,ls_activity,conf);
            hold off;
        end
    end

    if strcmp('png',p.mode)
        print([p.outfile,'.png'],'-dpng','-S640,480');
    elseif strcmp('paper',p.mode)
        print([p.outfile,'.eps'],'-deps','-S320,240');
    elseif ~strcmp('monitor',p.mode)
        error('%s: unkown print mode: %s!',upper(mfilename),p.mode);
    end

else

    % === Plot the wave field using Gnuplot ===
    % Data file name
    datafile = 'wavefield.dat';

    % Save the data for plotting with Gnuplot
    % NOTE: we have to transpose the wave field P, because the x and y values
    % are stored in a different way for Gnupot.
    gp_save_matrix(datafile,x,y,real(P'))

    % Check if we have used a loudspeaker array that is larger than the listener
    % area. If so, reduce the array length in order to avoid Gnuplot from
    % drawing loudspeakers out of the graph.
    % FIXME: check this for the different array forms we have!
    if L>abs(x(1)-x(end))
        L = abs(x(1)-x(end));
    end
    % Loudspeaker positions and directions
    [x0,y0,phi] = secondary_source_positions(L,conf);
    nLS = length(x0);

    % Check if we should plot the loudspeakers. If not set their distance to a
    % value <0.1 in order to avoid the plotting by Gnuplot
    if(p.loudspeakers)
        gp_LSdist = LSdist;
        if  LSdist<= 0.01
            warning(['%s: the given loudspeaker distance is to small. ',...
                    'Disabling plotting of the loudspeakers'],upper(mfilename));
        end
    else
        gp_LSdist = 0.01;
    end

    if strcmp('paper',p.mode)
    % Generate the Gnuplot command line
        cmd = sprintf(['gnuplot<<EOC\n', ...
            'set loadpath "%s/SFS_plotting"\n', ...
            'call "gplatex_plot_wavefield.gnu" "%f" "%f" "%f" "%f" ', ...
        '"%f" "%f" "%f" "%s" "%s"\nEOC\n', ...
            'epstopdf %s.eps'], ...
            sfspath,x(1),x(end),y(1),y(end),nLS,gp_LSdist,x0(1),datafile, ...
        outfile,outfile);
    elseif strcmp('monitor',p.mode)
        cmd = sprintf(['gnuplot<<EOC\n', ...
            'set loadpath "%s/SFS_plotting"\n', ...
            'call "gp_plot_wavefield.gnu" "%f" "%f" "%f" "%f" ', ...
        '"%f" "%f" "%f" "%s"\nEOC\n'], ...
        sfspath,x(1),x(end),y(1),y(end),nLS,gp_LSdist,x0(1),datafile);
    elseif strcmp('png',p.mode)
        to_be_implemented;
    else
        error('%s: %s is not a valid plotting mode!',upper(mfilename),p.mode);
    end

    % Start Gnuplot for plotting the data
    system(cmd);

    % Remove data file
    delete(datafile);

end
