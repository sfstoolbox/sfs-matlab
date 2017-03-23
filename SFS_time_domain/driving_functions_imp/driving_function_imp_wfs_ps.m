function [delay,weight] = driving_function_imp_wfs_ps(x0,nx0,xs,conf)
%DRIVING_FUNCTION_IMP_WFS_PS calculates the WFS weighting and delaying for a
%point source as source model
%
%   Usage: [delay,weight] = driving_function_imp_wfs_ps(x0,nx0,xs,conf)
%
%   Input parameters:
%       x0      - position  of secondary sources / m [nx3]
%       nx0     - direction of secondary sources [nx3]
%       xs      - position of point source [nx3]
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       delay   - delay of the driving function / s
%       weight  - weight (amplitude) of the driving function
%
%   DRIVING_FUNCTION_IMP_WFS_PS(x0,nx0,xs,conf) returns delays and weights for
%   the WFS driving function for a point source as source model.
%
%   See also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs_ps

%   References:
%       E. Start (1997) - "Direct Sound Enhancement by Wave Field Synthesis", 
%       PhD thesis, TU Delft
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,xs);
isargstruct(conf);


%% ===== Configuration ==================================================
% Speed of sound
c = conf.c;
xref = conf.xref;
fs = conf.fs;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation =====================================================
% Get the delay and weighting factors
if strcmp('2D',dimension) || strcmp('3D',dimension)

    % === 2- or 3-Dimensional ============================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % d using a point source as source model
        %
        %                   1  (x0-xs) nx0
        % d(x0,t) = h(t) * --- ----------- delta(t-|x0-xs|/c)
        %                  2pi  |x0-xs|^2
        %
        % See http://sfstoolbox.org/#equation-d.wfs.ps
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Delay and amplitude weight
        delay = 1./c .* r;
        weight = 1./(2.*pi) .* vector_product(x0-xs,nx0,2) ./ r.^2;
        %
    case 'legacy'
        % --- Old SFS Toolbox default ------------------------------------
        % d using a point source as source model
        %
        %                   1   (x0-xs) nx0
        % d(x0,t) = h(t) * --- ------------- delta(t-|x0-xs|/c)
        %                  2pi |x0-xs|^(3/2)
        %
        % See Wierstorf (2014), eq.(2.63)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Delay and amplitude weight
        delay = 1/c .* r;
        weight = 1/(2*pi) .* vector_product(x0-xs,nx0,2) ./ r.^(3/2);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

     % === 2.5-Dimensional ================================================

     % Reference point
     xref = repmat(xref,[size(x0,1) 1]);
     switch driving_functions
     case {'default', 'reference_point'}
         % Driving function with only one stationary phase approximation,
         % reference to one point in field
         %
         % r = |x0-xs|
         r = vector_norm(x0-xs,2);
         % 2.5D correction factor
         %         _____________________
         %        |      |xref-x0|
         % g0 = _ |---------------------
         %       \| |xref-x0| + |x0-xs|
         %
         % See Start (1997), eq. (3.11)
         %
         g0 = sqrt( vector_norm(xref-x0,2) ./ (vector_norm(x0-xref,2) + r) );
         %                                 ___
         %                                | 1    (x0-xs) nx0
         % d_2.5D(x0,t) = h_pre(t) * g0 _ |---  ------------- delta(t-|x0-xs|/c)
         %                               \|2pi  |x0-xs|^(3/2)
         %
         % See http://sfstoolbox.org/#equation-d.wfs.ps.2.5D
         %
         % Delay and amplitude weight
         delay = 1./c .* r;
         weight = g0 ./ sqrt(2.*pi) .* vector_product(x0-xs,nx0,2) ./ r.^(3/2);
         %
     case 'reference_line'
         % Driving function with two stationary phase approximations,
         % reference to a line parallel to a LINEAR secondary source distribution
         %
         % Distance ref-line to linear ssd
         dref = vector_product(xref-x0,nx0,2);
         % Distance source and linear ssd
         ds = vector_product(xs-x0,nx0,2);
         %
         % 2.5D correction factor
         %        _______________________
         % g0 = \| d_ref / (d_ref - d_s)
         %
         % See Start (1997), eq. (3.16)
         %
         g0 = sqrt( dref ./ (dref - ds) );
         %                                 ___
         %                                | 1    (x0-xs) nx0
         % d_2.5D(x0,t) = h_pre(t) * g0 _ |---  ------------- delta(t-|x0-xs|/c)
         %                               \|2pi  |x0-xs|^(3/2)
         %
         % Inverse Fourier Transform of Start (1997), eq. (3.17)
         %
         % r = |x0-xs|
         r = vector_norm(x0-xs,2);
         % Delay and amplitude weight
         delay = 1./c .* r;
         weight = g0 ./ sqrt(2.*pi) .* vector_product(x0-xs,nx0,2) ./ r.^(3./2);
         %
     case 'legacy'
         % --- SFS Toolbox ------------------------------------------------
         % 2.5D correction factor
         %        ______________
         % g0 = \| 2pi |xref-x0|
         %
         g0 = sqrt(2*pi*vector_norm(xref-x0,2));
         %
         % d_2.5D using a point source as source model
         %
         %                        g0  (x0-xs) nx0
         % d_2.5D(x0,t) = h(t) * --- ------------- delta(t-|x0-xs|/c)
         %                       2pi |x0-xs|^(3/2)
         %
         % See Wierstorf (2014), eq.(2.64)
         %
         % r = |x0-xs|
         r = vector_norm(x0-xs,2);
         % Delay and amplitude weight
         delay = 1./c .* r;
         weight = g0 ./ (2.*pi) .* vector_product(x0-xs,nx0,2) ./ r.^(3./2);
         %
     otherwise
         error(['%s: %s, this type of driving function is not implemented', ...
           'for a 2.5D point source.'],upper(mfilename),driving_functions);
     end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
