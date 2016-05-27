function D = driving_function_mono_wfs_mps(x0,nx0,xs,vs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_PS returns the driving signal D for a point source in
%WFS
%
%   Usage: D = driving_function_mono_wfs_mps(x0,nx0,xs,vs,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       xs          - current position of virtual moving point source / m [nx3]
%       vs          - velocity vector of moving point source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_MPS(x0,nx0,xs,vs,f,conf) returns WFS driving
%   signals for the given secondary sources, the virtual moving point source
%   position and the frequency f.
%
%   References:
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
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
nargmin = 6;
nargmax = 6;
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

% Frequency
omega = 2*pi*f;


if strcmp('2D',dimension) || strcmp('3D',dimension)
    
    % === 2- or 3-Dimensional ============================================
    
    
    
elseif strcmp('2.5D',dimension)
    
    % === 2.5-Dimensional ================================================
    
    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    switch driving_functions
        case {'default', 'reference_point'}
            % Source model for a uniformly moving point source: retarded 3D Green's
            % function.
            %
            %                 1   e^(+-i w/c R)
            % G(x-xs,vs,w) = --- --------------
            %                4pi      R_1
            %
            % vs = v * ns   % velocity vector
            % M = v/c  % Mach number
            % xparallel = (x-xs)*ns  % component in direction of movement
            % R_1 = sqrt( M^2*xparallel.^2 + (1-M^2)*|x-xs|.^2 )
            % R = (M*xparallel + R1)./(1-M*2);
            %
            
            [tau, R, R1] = retarded_time(x0-xs, 0, vs, conf);
  
            % 2.5D correction factor
            %         _____________________
            %        |      |xref-x0|
            % g0 = _ |---------------------
            %       \| |xref-x0| +  R (x
            %
            g0 = sqrt( vector_norm(xref-x0,2)./ (vector_norm(x0-xref,2) + R));
            %
            % D_2.5D(x0,w) =
            %       ___    ___
            %      | 1    |i w (x0-xs) nx0
            % g0 _ |--- _ |--- ------------- e^(-i w0/c |x0-xs|)
            %     \|2pi  \| c    R1^(3/2)
            %
            % Driving signal
            D = sqrt(1i*omega/(2*pi*c)) .* g0 .* vector_product(x0-xs,nx0,2) ...
                ./ R1.^(3/2) .* exp(-1i*omega*tau);
            %
        otherwise
            error(['%s: %s, this type of driving function is not ', ...
                ' implemented for a 2.5D point source.'], upper(mfilename), ...
                driving_functions);
    end
    
else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
