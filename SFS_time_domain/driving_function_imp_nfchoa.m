function [d] = driving_function_imp_nfchoa(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_NFCHOA calculates the NFC-HOA driving function
%
%   Usage: [d] = driving_function_imp_nfchoa(x0,xs,src,[conf]);
%
%   Input parameters:
%       x0      - position  and direction of secondary sources / m
%       xs      - position of virtual source or direction of plane wave / m
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       d  - matrix of driving signals
%
%   DRIVING_FUNCTION_IMP_NFCHOA(x0,xs,src,conf) returns the
%   driving function of NFC-HOA for the given source type and position,
%   and loudspeaker positions.
%
%   see also: driving_function_imp_nfchoa_ps, sound_field_imp_nfchoa

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
isargsecondarysource(x0)
isargxs(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
nls = size(x0,1);
N = conf.N;
X0 = conf.secondary_sources.center;


%% ===== Computation =====================================================

% generate stimulus pusle
pulse = dirac_imp();
% radius of array
R = norm(x0(1,1:3)-X0);
% get maximum order of spherical harmonics
order = nfchoa_order(nls,conf);

% if-request as a workaround for the right direction of the sound field
if strcmpi(src,'pw')
    [theta_src, r_src] = cart2pol(-xs(1),xs(2));
elseif strcmpi(src,'ps')
    [theta_src, r_src] = cart2pol(xs(1),-xs(2));
else
    [theta_src, r_src] = cart2pol(xs(1),xs(2));
end

% compute impulse responses of modal filters
dm = zeros(order+1,N);
for n=1:order+1
    dm(n,:) = [pulse zeros(1,N-length(pulse))];
    
    % get the second-order sections for the different virtual sources
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

    % apply them by a bilinear transform and filtering
    [b,a] = bilinear_transform(sos,conf);
    for ii=1:length(b)
        dm(n,:) = filter(b{ii},a{ii},dm(n,:));
    end
end

% compute input signal for IFFT
d = zeros(2*order+1,N);
for n=-order:order
    d(n+order+1,:) = dm(abs(n)+1,:) .* exp(1i*n*theta_src);
end

if(iseven(nls))
   d = d(2:end,:);
end

% spatial IFFT
d = circshift(d,[order+1 0]);
d = (2*order+1)*ifft(d,[],1);
d = real(d');

% subsample d if we have fewer secondary sources than the applied order
if size(d,2)>nls
    % check if we have a multiple of the order
    if mod(size(d,2),nls)~=0
        conf_tmp = conf;
        conf_tmp.nfchoa.order = [];
        error(['%s: the given number of driving signals (%i) can not ', ...
            'be subsampled to %i secondary sources. Choose a NFC-HOA ', ...
            'order that is a multiple of %i.'], ...
            upper(mfilename),size(d,2),nls,nfchoa_order(nls,conf_tmp));
    end
    % subsample d
    ratio = size(d,2)/nls;
    d = d(:,1:ratio:end);
% subsample the secondary sources if we have fewer driving signals than
% secondary sources
elseif size(d,2)<nls
    % check if we have a multiple of the secondary sources
    if mod(nls,size(d,2))~=0
        conf_tmp = conf;
        conf_tmp.nfchoa.order = [];
         error(['%s: the given number of secondary sources (%i) can not ', ...
            'be subsampled to %i driving signals. Choose a NFC-HOA ', ...
            'order that is a multiple of %i.'], ...
            upper(mfilename),nls,size(d,2),nfchoa_order(size(d,2)));
    end
    % subsample x0
    ratio = nls/size(d,2);
    d_new = zeros(N,nls);
    d_new(:,1:ratio:end) = d;
    d = d_new;
end
