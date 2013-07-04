function cmd = gp_get_plot_command(p)
%GP_GET_PLOT_COMMAND returns the gnuplot command for plotting the sound field
%and loudspeaker
%
%   Usage: cmd = gp_get_plot_command(p)
%
%   Input parameters:
%       p   - struct holding all the conf.plot data and more
%
%   Output parameters:
%       cmd - plot command to insert in gnuplot file
%
%   GP_GET_PLOT_COMMAND(p) returns the plot command for gnuplot to plot the
%   sound field and loudspeaker as symbols or points, or completly without
%   loudspeakers.
%
%   See also: gp_print_*, plot_wavefield

if p.loudspeakers && p.realloudspeakers
    % Plotting real loudspeaker symbols
    cmd = sprintf([ ...
        'call ''gp_draw_loudspeakers.gnu'' ''%s'' ''%f''\n', ...
        'plot ''%s'' binary matrix with image'], ...
        p.lsfile,p.lssize,p.datafile);
elseif p.loudspeakers
    % Plotting only points at the loudspeaker positions
    cmd = sprintf([ ...
        'plot ''%s'' binary matrix with image,\\\n', ...
        '     ''%s'' u 1:2 w points ls 1'], ...
        p.datafile,p.lsfile);
else
    % Plotting no loudspeakers at all
    cmd = sprintf( ...
        'plot ''%s'' binary matrix with image', ...
        p.datafile);
end
