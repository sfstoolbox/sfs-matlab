function [delay,weight] = driving_function_imp_wfs_pw(x0,nx0,nk,conf)
%DRIVING_FUNCTION_IMP_WFS_PW calculates the WFS weighting and delaying for a
%plane wave as source model
%
%   Usage: [delay,weight] = driving_function_imp_wfs_pw(x0,nx0,nk,[conf]);
%
%   Input parameters:
%       x0      - position  of secondary sources (m) [nx3]
%       nx0     - direction of secondary sources [nx3]
%       nk      - direction of plane wave [nx3]
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       delay   - delay of the driving function (s)
%       weight  - weight (amplitude) of the driving function
%
%   DRIVING_FUNCTION_IMP_WFS_PW(x0,nx0,nk,conf) returns delays and weights for
%   the WFS driving function for plane wave as source model.
%   
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, Tu Berlin
%
%   see also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs_pw

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,nk);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


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

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % d_2D using a plane wave as source model
        %
        % d_2D(x0,t) = h(t) * -2 nk nx0 delta(t - 1/c nk x0)
        %
        % see Wierstorf (2014), p.25 (2.45)
        %
        % Delay and amplitude weight
        delay = 1/c * vector_product(nk,x0,2);
        weight = -2 .* vector_product(nk,nx0,2);
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % Reference point
        xref = repmat(xref,[size(x0,1) 1]);
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % d_2.5D using a plane wave as source model
        %
        % d_2.5D(x0,t) = h(t) * -2 g0 nk nx0 delta(t - 1/c nk x0)
        % 
        % see Wierstorf (2014), p.25 (2.46)
        %
        % Delay and amplitude weight
        delay = 1/c .* vector_product(nk,x0,2);
        weight = -2*g0 .* vector_product(nk,nx0,2);
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
