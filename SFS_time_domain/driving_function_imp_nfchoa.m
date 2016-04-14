function [d,dm] = driving_function_imp_nfchoa(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_NFCHOA calculates the NFC-HOA driving function
%
%   Usage: [d] = driving_function_imp_nfchoa(x0,xs,src,conf)
%
%   Input parameters:
%       x0      - position  and direction of secondary sources / m
%       xs      - position of virtual source or direction of plane wave / m
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       d  - matrix of driving signals
%       dm - matrix of driving funtion in spherical/circular domain
%
%   DRIVING_FUNCTION_IMP_NFCHOA(x0,xs,src,conf) returns the
%   driving function of NFC-HOA for the given source type and position,
%   and loudspeaker positions.
%
%   See also: driving_function_imp_nfchoa_ps, sound_field_imp_nfchoa

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
nls = size(x0,1);
N = conf.N;
X0 = conf.secondary_sources.center;

%% ===== Computation =====================================================
% Generate stimulus pusle
pulse = dirac_imp();
% Radius of array
R = norm(x0(1,1:3)-X0); 
% Ambisonics order
if isempty(conf.nfchoa.order)
    % Get maximum order of spherical harmonics
    order = nfchoa_order(nls,conf);
else
    order = conf.nfchoa.order;
end

% Correct position of source for off-center arrays
xs(1:3) = xs(1:3)-X0;
[theta_src, r_src] = cart2pol(xs(1),xs(2));

% Compute impulse responses of modal filters
dm = [repmat(pulse,[order+1 1]) zeros(order+1,N-length(pulse))];
for n=2:order+1

    % Get the second-order sections for the different virtual sources
    if strcmp('pw',src)
        % === Plane wave =================================================
        sos = driving_function_imp_nfchoa_pw(n-1,R,conf);
    elseif strcmp('ps',src)
        % === Point source ===============================================
        sos = driving_function_imp_nfchoa_ps(n-1,R,r_src,conf);
    elseif strcmp('ls',src)
        % === Line source ================================================
        sos = driving_function_imp_nfchoa_ls(n-1,R,r_src,conf);
    elseif strcmp('fs',src)
        % === Focussed source ============================================
        sos = driving_function_imp_nfchoa_fs(n-1,R,r_src,conf);
    else
        error('%s: %s is not a known source type.',upper(mfilename),src);
    end

    % Apply them by a bilinear transform and filtering
    [b,a] = bilinear_transform(sos,conf);
    for ii=1:length(b)
        dm(n,:) = filter(b{ii},a{ii},dm(n,:));
    end
end

% -------------------------------------------------------------------------
%                      ___
%                1     \
% D(phi0,w) = -------  /__     Hm(w) e^(im(phi0-phi_src))
%             2 pi r0  m=-M..M
%
% see Spors et al. (2011), eq.(4) and (5)
%--------------------------------------------------------------------------

% Compute input signal for IFFT
dM = zeros(2*order+1,N);
for n=-order:order
    dM(n+order+1,:) = dm(abs(n)+1,:) * exp(-1i*n*theta_src);
end
% spatial IFFT
d = zeros(nls,N);
for l=1:N
    d(:,l) = sum(buffer(dM(:,l),nls),2);
end
d = circshift(d,[-order,0]);
d = ifft(transpose(d),[],2);
d = 1/pi/R*nls*real(d);

% -------------------------------------------------------------------------
% The following is the direct implementation of the spatial IDFT which
% takes longer time for higher orders.
% 
if(0)
    d = zeros(N,nls);
    for n=1:nls
        phin = 2*pi/nls*(n-1); % first phi0 always 0 ??
        dtemp = zeros(1,N);
        for m=-order:order
            dtemp = dtemp + dm(abs(m)+1,:) .* exp(1i*m*(phin-theta_src));
        end
        d(:,n) = transpose(dtemp);
    end
    d = 1/pi/R*real(d);
end