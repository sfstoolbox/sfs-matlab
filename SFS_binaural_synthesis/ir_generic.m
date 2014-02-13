function ir = ir_generic(X,phi,x0,d,irs,conf)
%IR_GENERIC Generate a IR
%
%   Usage: ir = ir_generic(X,phi,x0,d,irs,[conf])
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%                 0 means the head is oriented towards the x-axis.
%       x0      - secondary sources [n x 6] / m
%       d       - driving signals [m x n]
%       irs     - IR data set for the secondary sources
%       conf    - optional configuration struct (see SFS_config) 
%
%   Output parameters:
%       ir      - Impulse response for the desired driving functions (nx2 matrix)
%
%   IR_GENERIC(X,phi,x0,d,irs,conf) calculates a binaural room impulse
%   response for the given secondary sources and driving signals.
%
%   see also: ir_wfs, ir_nfchoa, ir_point_source, auralize_ir

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
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(X);
isargscalar(phi);
isargsecondarysource(x0);
isargmatrix(d);
check_irs(irs);
if nargin==nargmax-1
    conf = SFS_config;
end


%% ===== Configuration ==================================================
N = conf.N;                   % target length of BRS impulse responses


%% ===== Variables ======================================================
phi = correct_azimuth(phi);


%% ===== BRIR ===========================================================
% Initial values
ir_generic = zeros(N,2);

% Create a BRIR for every single loudspeaker
warning('off','SFS:irs_intpol');
for ii=1:size(x0,1)

    % direction of the source from the listener
    x_direction = x0(ii,1:3)-X;
    % change to spherical coordinates
    [alpha,theta,r] = cart2sph(x_direction(1),x_direction(2),x_direction(3));

    % === Secondary source model: Greens function ===
    g = 1./(4*pi*r);

    % Incoporate head orientation and ensure -pi <= alpha < pi
    alpha = correct_azimuth(alpha-phi);

    % === IR interpolation ===
    % Get the desired IR.
    % If needed interpolate the given IR set
    ir = get_ir(irs,[alpha,theta,r],conf);

    % === Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay ===
    % Also applying the weights of the secondary sources including integration
    % weights or tapering windows etc.
    ir_generic = ir_generic + fix_length(convolution(ir,d(:,ii)),N).*g.*x0(ii,7);

end
warning('on','SFS:irs_intpol');


%% ===== Headphone compensation =========================================
ir = compensate_headphone(ir_generic,conf);
