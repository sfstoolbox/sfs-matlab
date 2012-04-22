% formats MATLAB figures for various publication types
% S.Spors / 23.06.04

% FIXME: is this function needed anymore?
function GraphDefaults(p_type)

% FIXME: this needs documentation for the avaiable figure sizes.
% TODO: but is a nice idea and should also be implemented with Gnuplot in order
% to save a lot of time!

switch lower(p_type)

    % two figures side by side on one page
    case 'two_diss'
        figsize(7,7);
        fontsize(10);
        fontname('Times');

    % two figures side by side on one page with colorbar
    case 'two_diss2'
        figsize(6.7,5.0);
        fontsize(10);
        fontname('Times');

    % one figure by side, full page
    case 'full_diss'
        figsize(16,22);
        fontsize(10);
        fontname('Times');

    % one figure by side, half page
    case 'half_diss'
        figsize(16,10);
        fontsize(10);
        fontname('Times');

    case 'half_diss2'
        figsize(16,8);
        fontsize(10);
        fontname('Times');

    % 8polar plots, half page
    case '8polar'
        figsize(26,10);
        fontsize(10);
        fontname('Times');

    % one figure by side, 1/3rd page height
    case 'one_diss'
        figsize(15,7);
        fontsize(10);
        fontname('Times');

    % figure with box shape, about 1/3rd page height
    case 'box_diss'
        figsize(10,10);
        fontsize(10);
        fontname('Times');

    case 'talk'
        figsize(25,14);
        fontsize(18);
        fontname('Times');

    case 'two_talk'
        figsize(12,12);
        fontsize(18);
        fontname('Times');

    case 'box_talk'
        set(0,'DefaultLineLineWidth',2);
        figsize(14,14);
        fontsize(18);
        fontname('Times');

    case 'movie16_9'
        set(0,'DefaultLineLineWidth',2);
        figsize(15.31,7.74);
        fontsize(16);
        fontname('Times');

     case 'movie'
        set(0,'DefaultLineLineWidth',2);
        figsize(10.98,7.74);
        fontsize(16);
        %fontname('Times');

    case 'paper'
        set(0,'DefaultLineLineWidth',2);
        figsize(8,6);
        %figsize(12,9);
        fontsize(10);
        fontname('Times');

    case 'paper2'
        figsize(10,5);
        fontsize(10);
        fontname('Times');

    case 'png'
        set(0,'DefaultLineLineWidth',2);
        figsize(12,9);
        fontsize(10);
        fontname('Helvetica');

    case 'monitor'
        set(0,'DefaultLineLineWidth',2);
        figsize(12,9);
        fontsize(10);
        fontname('Helvetica');

    case 'hagen'
        set (0,'defaultlinecolor',[0 0.37647 0.67843]);

    otherwise
        error('%s: unknown figure format!',upper(mfilename));

end
