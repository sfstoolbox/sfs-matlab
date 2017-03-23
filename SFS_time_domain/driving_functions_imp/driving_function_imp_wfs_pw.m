function [delay,weight] = driving_function_imp_wfs_pw(x0,nx0,nk,conf)
%DRIVING_FUNCTION_IMP_WFS_PW calculates the WFS weighting and delaying for a
%plane wave as source model
%
%   Usage: [delay,weight] = driving_function_imp_wfs_pw(x0,nx0,nk,conf);
%
%   Input parameters:
%       x0      - position  of secondary sources / m [nx3]
%       nx0     - direction of secondary sources [nx3]
%       nk      - direction of plane wave [nx3]
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       delay   - delay of the driving function / s
%       weight  - weight (amplitude) of the driving function
%
%   DRIVING_FUNCTION_IMP_WFS_PW(x0,nx0,nk,conf) returns delays and weights for
%   the WFS driving function for plane wave as source model.
%
%   See also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs_pw

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
isargmatrix(x0,nx0,nk);
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
        % d_2D using a plane wave as source model
        %
        % d_2D(x0,t) = h(t) * 2 nk nx0 delta(t - 1/c nk x0)
        %
        % See http://sfstoolbox.org/#equation-d.wfs.pw
        %
        % Delay and amplitude weight
        delay = 1./c .* vector_product(nk,x0,2);
        weight = 2 .* vector_product(nk,nx0,2);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    switch driving_functions
    case {'default', 'reference_point'}      
        % Driving function with only one stationary phase approximation, i.e.
        % reference to one point in field
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt( 2*pi*vector_norm(xref-x0,2) );
        %
        % d_2.5D using a plane wave as source model
        %
        % d_2.5D(x0,t) = h(t) * 2 g0 nk nx0 delta(t - 1/c nk x0)
        %
        % See http://sfstoolbox.org/en/update_wfs_ps/#equation-d.wfs.pw.2.5D
        %
        % Delay and amplitude weight
        delay = 1./c .* vector_product(nk,x0,2);
        weight = 2.*g0 .* vector_product(nk,nx0,2);
        %
    case 'reference_line'
        % Driving function with two stationary phase approximations,
        % reference to a line parallel to a LINEAR secondary source distribution
        %
        % Distance ref-line to linear ssd
        dref = vector_product(xref-x0,nx0,2);
        %
        % 2.5D correction factor
        %        ____________
        % g0 = \| 2pi * d_ref
        %
        g0 = sqrt(2.*pi.*dref);
        %
        % d_2.5D using a plane wave as source model
        %
        %                             ______   
        % d_2.5D(x0,w) = h(t) * 2g0 \|nk nx0 delta(t - 1/c nk x0)
        % 
        %
        % See Schultz (2016), eq. (2.183)
        %
        % Delay and amplitude weight
        delay = 1./c .* vector_product(nk,x0,2);
        weight = 2.*g0 .* sqrt(vector_product(nk,nx0,2));
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
