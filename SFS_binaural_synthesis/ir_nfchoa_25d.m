function ir = ir_nfchoa_25d(X,phi,xs,src,L,irs,conf)
%IR_NFCHOA_25D Generate a IR for HOA
%
%   Usage: ir = ir_nfchoa_25d(X,phi,xs,src,L,irs,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%                 0 means the head is oriented towards the x-axis.
%       xs      - virtual source position [ys > Y0 => focused source] (m)
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%       L       - Length of loudspeaker array (m)
%       irs     - IR data set for the secondary sources
%       conf    - optional configuration struct (see SFS_config) 
%
%   Output parameters:
%       ir      - Impulse response for the desired HOA synthesis (nx2 matrix)
%
%   IR_NFCHOA_25D(X,phi,xs,src,L,irs,conf) calculates a binaural room impulse
%   response for a virtual source at xs for a virtual HOA array and a
%   listener located at X.
%
%   see also: brs_nfchoa_25d, ir_wfs_25d, ir_point_source, auralize_ir

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
[X,xs] = position_vector(X,xs);
if conf.debug
    isargscalar(phi);
    isargpositivescalar(L);
    isargchar(src);
    check_irs(irs);
end


%% ===== Configuration ==================================================
N = conf.N;                   % target length of BRS impulse responses


%% ===== Variables ======================================================
phi = correct_azimuth(phi);
% Loudspeaker positions
x0 = secondary_source_positions(L,conf);
nls = size(x0,1);


%% ===== BRIR ===========================================================
% Initial values
ir_hoa = zeros(N,2);

% calculate driving function
d = driving_function_imp_nfchoa_25d(x0,xs,src,L,conf);

% Create a BRIR for every single loudspeaker
warning('off','SFS:irs_intpol');
for ii=1:nls

    % === Secondary source model: Greens function ===
    g = 1./(4*pi*norm(X-x0(ii,1:3)));

    % === Secondary source angle ===
    % Calculate the angle between the given loudspeaker and the listener.
    % This is needed for the HRIR dataset.
    %
    %                             y-axis
    %                               ^
    %                               |
    %                               |
    %                               |
    %            [X Y], phi=0       |
    %              O------------    |  a = alpha
    %               \ a |           |  tan(alpha) = (y0-Y)/(x0-X)
    %                \ /            |
    %                 \             |
    %                  \            |
    %   -------v--v--v--v--v--v--v--v--v--v--v--v--v--v--v------> x-axis
    %                [x0 y0]
    %
    % Angle between listener and secondary source (-pi < alpha <= pi)
    % Note: phi is the orientation of the listener (see first graph)
    [alpha,theta_tmp,r_tmp] = cart2sph(x0(ii,1)-X(1),x0(ii,2)-X(2),0);
    %
    % Ensure -pi <= alpha < pi
    alpha = correct_azimuth(alpha-phi);

    % === IR interpolation ===
    % Get the desired IR.
    % If needed interpolate the given IR set
    ir = get_ir(irs,alpha,0);

    % === Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay ===
    ir_hoa(:,1) = ir_hoa(:,1) + fix_ir_length(conv(ir(:,1),d(:,ii)),N) .* g;
    ir_hoa(:,2) = ir_hoa(:,2) + fix_ir_length(conv(ir(:,2),d(:,ii)),N) .* g;

end
warning('on','SFS:irs_intpol');


%% ===== Headphone compensation =========================================
ir = compensate_headphone(ir_hoa,conf);
