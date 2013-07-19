function D = driving_function_mono_wfs_ls(x0,nx0,xs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_LS returns the driving signal D for a line source in
%WFS
%
%   Usage: D = driving_function_mono_wfs_ls(x0,nx0,xs,f,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       xs          - position of virtual line source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_LS(x0,xs,f,src,conf) returns WFS driving signals
%   for the given secondary sources, the virtual line source position and the
%   frequency f.
%
%   References:
%       FIXME: Update references
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: driving_function_mono_wfs, driving_function_imp_wfs_ps

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
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% frequency
omega = 2*pi*f;


if strcmp('2D',dimension)
    
    % === 2-Dimensional ==================================================
    
    % Ensure 2D
    x0 = x0(:,1:2);
    nx0 = nx0(:,1:2);
    xs = xs(:,1:2);
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % D_2D using a line source
        %
        %                 iw (x0-xs) nx0   (2)/ w         \
        % D_2D(x0,w) =  - -- -----------  H1  | - |x0-xs| |
        %                 2c   |x0-xs|        \ c         /
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1i*omega/(2*c) .* vector_product(x0-xs,nx0,2) ./ r .* ...
            besselh(1,2,omega/c.*r);
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D line source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)
    
    % === 2.5-Dimensional ================================================
    
    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a line source
        %                        ___
        %                   1   |i w  (x0-xs) nx0   (2)/ w         \
        % D_2.5D(x0,w) =  - - _ |---  -----------  H1  | - |x0-xs| |
        %                   2  \| c    |x0-xs|         \ c         /
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1/2.*g0 .* sqrt(i*omega/c) .* vector_product(x0-xs,nx0,2) ./ r .* ...
            besselh(1,2,omega/c.*r);
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D line source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)
    

    % === 3-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        warning(['%s: you use conf.dimension="3D" together with a line ', ...
            'which will give no meaningfull results.',upper(mfilename));
        % --- SFS Toolbox ------------------------------------------------
        % D_3D using a line source
        %
        %                 iw (x0-xs) nx0   (2)/ w         \
        % D_3D(x0,w) =  - -- -----------  H1  | - |x0-xs| |
        %                 2c   |x0-xs|        \ c         /
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1i*omega/(2*c) .* vector_product(x0-xs,nx0,2) ./ r .* ...
            besselh(1,2,omega/c.*r);
        %
        to_be
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D line source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
