function D = driving_function_mono_wfs_ls(x0,nx0,xs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_LS returns the driving signal D for a line source in
%WFS
%
%   Usage: D = driving_function_mono_wfs_ls(x0,nx0,xs,nxs,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       xs          - position and orientation of virtual line source / m [nx3]
%                     or [nx6]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_LS(x0,nx0,xs,f,src,conf) returns WFS driving
%   signals for the given secondary sources, the virtual line source position,
%   its orientation xs(:,4:6), which is parallel to the line source, and the
%   frequency f. If no explicit orientation is given, [0 0 1] is assumed.
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,xs);
isargpositivescalar(f);
isargstruct(conf);


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

omega = 2*pi*f;
[xs,nxs] = get_position_and_orientation_ls(xs,conf);

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % D using a line source
        %
        %              iw (x0-xs) nx0   (2)/ w         \
        % D(x0,w) =  - -- -----------  H1  | - |x0-xs| |
        %              2c   |x0-xs|        \ c         /
        %
        % See http://sfstoolbox.org/#equation-D.wfs.ls
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -1i.*omega./(2.*c) ...
            .* vector_product(x0-xs,nx0,2) ./ r ...
            .* besselh(1,2,omega./c.*r);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a line source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a line source
        %                         ___
        %                   g0   |i w  (x0-xs) nx0   (2)/ w         \
        % D_2.5D(x0,w) =  - -- _ |---  -----------  H1  | - |x0-xs| |
        %                   2   \| c    |x0-xs|         \ c         /
        %
        % See http://sfstoolbox.org/#equation-D.wfs.ls.2.5D
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -g0./2 .* sqrt(i.*omega./c) ...
            .* vector_product(x0-xs,nx0,2) ./ r ...
            .* besselh(1,2,omega./c.*r);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D line source.'],upper(mfilename),driving_functions);
    end

elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % D using a line source
        %
        %              iw   v nx0    (2)/ w     \
        % D(x0,w) =  - -- --------  H1  | - |v| | ,
        %              2c   |v|         \ c     /
        %
        % where v = x0-xs - <x0-xs,nxs > nxs,
        % and |nxs| = 1.
        %
        % See http://sfstoolbox.org/#equation-d.wfs.ls
        % and http://sfstoolbox.org/#equation-v.ls
        %
        % v = (I - nxs'nxs)(x0-xs)
        % r = |v|
        nxs = nxs(1,:);
        v = (x0 - xs)*(eye(3) - nxs'*nxs);
        r = vector_norm(v,2);
        % Driving signal
        D = -1i*omega/(2*c) .* vector_product(v,nx0,2) ./ r .* ...
            besselh(1,2,omega/c.*r);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a line source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
