function [sos,g] = driving_function_imp_nfchoa_ps(N,R,r,conf)
%DRIVING_FUNCTION_IMP_NFCHOA_PS calculates the second-order section
%representation for a virtual point source in NFC-HOA
%
%   Usage: sos = driving_function_imp_nfchoa_ps(N,R,r,conf)
%
%   Input parameters:
%       N       - order of spherical hankel function
%       R       - radius of secondary source array / m
%       r       - distance of point source from array center / m
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       sos     - second-order section representation
%       g       - scalar gain factor
%
%   DRIVING_FUNCTION_IMP_NFCHOA_PS(N,R,r,conf) returns the second-order section
%   representation for the NFC-HOA driving function for a virtual point source
%   as source model.
%
%   References:
%       S. Spors, V. Kuscher, J. Ahrens (2011) - "Efficient realization of
%       model-based rendering for 2.5-dimensional near-field compensated higher
%       order Ambisonics", WASPAA, p. 61-64
%
%   See also: sound_field_imp, sound_field_imp_nfchoa, driving_function_imp_nfchoa

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
isargpositivescalar(N,R,r);
isargstruct(conf);


%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation =====================================================
% Find spherical hankel function zeros
z = sphbesselh_zeros(N);

% Get the delay and weighting factors
if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension) || strcmp('3D',dimension)

    % === 2.5- & 3-Dimensional ==========================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D using a point source as source model
        %
        [sos,g] = zp2sos(z*c/r,z*c/R,1,'up','none');
        g = g * R/r;
        %
        % Compare Spors et al. (2011), eq. (11)
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end
else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
