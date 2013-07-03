function irs_pw = extrapolate_farfield_hrtfset_3d(irs,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET_3D far-field extrapolation of a given 3D HRTF 
%                                dataset
%
%   Usage: irs_pw = extrapolate_farfield_hrtfset(irs,[conf])
%
%   Input parameters:
%       irs     - IR data set for the virtual secondary sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       irs_pw  - IR data set extra polated to conation plane wave IRs
%
%   EXTRAPOLATE_FARFIELD_HRTFSET_3D(IRS) generates a far-field extrapolated 
%   set of impulse responses, using the given irs set. Far-field means that 
%   the resulting impulse responses are plane waves. The extrapolation is 
%   done via 3D WFS.
%
%   see also: ir_point_source, get_ir, driving_function_imp_wfs_3d

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut für Nachrichtentechnik                    *
%                         Universität Rostock                                *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
% check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
else
%     isargstruct(conf);
end


%% ===== Configuration ===================================================
fs = conf.fs;                   % sampling frequency


%% ===== Variables ======================================================
conf.array = 'spherical';
conf.usetapwin = 0;
conf.wfs.usehpre = 1;
conf.usefracdelay = 0;
conf.xref = [0 0 0];
x0 = zeros(length(irs.source_position),8);
x0(:,1:3) = irs.source_position.';
x0(:,4:6) = direction_vector(x0(:,1:3),repmat(conf.xref,length(irs.source_position),1));
x0(:,7) = (irs.distance.^2)' .* cos(irs.apparent_elevation)';
% fix weight for northpole
x0(1,7) = x0(2,7);
x0(:,8) = weights_for_points_on_a_sphere_rectangle(irs.apparent_azimuth,...
          irs.apparent_elevation,irs.distance)';
conf.wfs.hprefhigh = aliasing_frequency(x0,conf);
conf.wfs.hpreflow = 1;
%% ===== Computation =====================================================
% get virtual secondary source positions
% x0_all = secondary_source_positions(L,conf);
x0_all = x0;
% Initialize new irs set
irs_pw = irs;
irs_pw.description = 'Extrapolated HRTF set containing plane waves';
irs_pw.left = zeros(size(irs_pw.left));
irs_pw.right = zeros(size(irs_pw.right));
irs_pw.distance = 'Inf';

% Generate a irs set for all given angles
for ii = 1:length(irs.apparent_azimuth)
disp([ii length(irs.apparent_azimuth)])
    % direction of plane wave
    xs = -irs.distance(ii).*[cos(irs.apparent_azimuth(ii)).*cos(irs.apparent_elevation(ii)) ...
           sin(irs.apparent_azimuth(ii)).*cos(irs.apparent_elevation(ii)) ...
           sin(irs.apparent_elevation(ii))];
       
    % calculate active virtual speakers
    x0 = secondary_source_selection(x0_all,xs,'pw');

    % generate tapering window
    win =ones(size(x0,1));

    % sum up contributions from individual virtual speakers
    for l=1:size(x0,1)
        % Driving function to get weighting and delaying
        [a,delay] = driving_function_imp_wfs_3d(x0(l,:),xs,'pw',conf);
        dt = delay*fs + irs.distance(l)/conf.c*fs;%round(irs.distance(l)/conf.c*fs);
        w=a*win(l);
        % get IR for the secondary source position
        [phi,theta,r] = cart2sph(x0(l,1),x0(l,2),x0(l,3));
        ir_tmp = get_ir(irs,phi,theta,r,conf.xref');
        % truncate IR length
        ir_tmp = fix_ir_length(ir_tmp,length(ir_tmp(:,1)),dt);
%       irr = fix_ir_length(ir_tmp(:,2),length(ir_tmp(:,2)),dt);
        irl = ir_tmp(:,1);
        irr = ir_tmp(:,2);
        % delay and weight HRTFs
        irs_pw.left(:,ii) = irs_pw.left(:,ii) + delayline(irl',dt,w,conf)';
        irs_pw.right(:,ii) = irs_pw.right(:,ii) + delayline(irr',dt,w,conf)';
    end

end

%% ===== Pre-equalization ===============================================
irs_pw.left = wfs_preequalization3d(irs_pw.left,conf);
irs_pw.right = wfs_preequalization3d(irs_pw.right,conf);
