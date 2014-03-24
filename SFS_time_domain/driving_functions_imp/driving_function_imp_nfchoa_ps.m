function sos = driving_function_imp_nfchoa_ps(N,R,r,conf)
%DRIVING_FUNCTION_IMP_NFCHOA_PS calculates the second-order section
%representation for a virtual point source in NFC-HOA
%
%   Usage: sos = driving_function_imp_nfchoa_ps(N,R,r,[conf]);
%
%   Input parameters:
%       N       - order of spherical hankel function
%       R       - radius of secondary source array / m
%       r       - distance of point source from array center / m
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       sos     - second-order section representation
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
%   see also: sound_field_imp, sound_field_imp_nfchoa, driving_function_imp_nfchoa

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
isargpositivescalar(N,R,r);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation =====================================================

% find spherical hankel function zeros
z = sphbesselh_zeros(N);

% Get the delay and weighting factors
if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D using a point source as source model
        %
        sos = zp2sos(z*c/r,z*c/R,1,'up','none');
        %
        % compare Spors et al. (2011)
        %
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 3D point source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
