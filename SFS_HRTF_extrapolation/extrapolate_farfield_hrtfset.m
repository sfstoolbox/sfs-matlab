function brs = extrapolate_farfield_hrtfset(phi,irs,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET far-field extrapolation of a given HRTF dataset
%
%   Usage: brs = extrapolate_farfield_hrtfset(phi,xs,L,src,irs,[conf])
%
%   Input parameters:
%       phi     - listener direction [head orientation] (rad)
%       irs     - IR data set for the virtual secondary sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set

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


% FIXME: return not a BRS but a irs in irs format

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargscalar(phi);
check_irs(irs);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
N = conf.N;                     % Target length of BRIR impulse responses
angles = rad(conf.brsangles);   % Angles for the BRIRs
fs = conf.fs;                   % sampling frequency

%% ===== Variables ======================================================

phi = correct_azimuth(phi);
R = irs.distance;
L = 2*R;
nls = length(irs.apparent_azimuth);
lenir = length(irs.left(:,1));

% get virtual loudspeaker positions from HRTF dataset
conf.array = 'circle';
conf.dx0=2*pi*R/nls;
conf.xref = [0 0 0];
conf.x0=zeros(nls,6);
conf.x0(:,1:2)=[R*cos(irs.apparent_azimuth) ; R*sin(irs.apparent_azimuth)]';
conf.x0(:,4:6)=direction_vector(conf.x0(:,1:3),repmat(conf.xref,nls,1));


x0 = secondary_source_positions(L,conf);

brs = zeros(N,2*length(angles));

Acorr=-1.7;                     % DAGA 2011 R=0.5m -> pw
Af=Acorr*sin(angles);


% append zeros or truncate IRs to target length
if(lenir<N)
    irs.left=cat(1,ir.left,zeros(N-lenir,2));
    irs.right=cat(1,ir.right,zeros(N-lenir,2));
else
    irs.left=irs.left(1:N,:);
    irs.right=irs.right(1:N,:);
end


%% ===== Computation =====================================================
% Generate a BRS set for all given angles
for ii = 1:length(angles)

    % variables
    brir = zeros(N,2);
    % FIXME: this works only for plane waves (xs/src is ignored)
    xs = -[cos(angles(ii)) sin(angles(ii))];

    % calculate active virtual speakers
    ls_activity = secondary_source_selection(x0,xs,'pw');

    % generate tapering window
    win = tapwin(L,ls_activity,conf);

    % sum up contributions from individual virtual speakers
    aidx=find(ls_activity>0);
    for l=aidx'
        % Driving function to get weighting and delaying
        [a,delay] = driving_function_imp_wfs_25d(x0(l,:),xs,'pw',conf);
        dt = delay*fs + round(R/conf.c*fs);
        w=a*win(l);

        % delay and weight HRTFs
        brir(:,1) = brir(:,1) + delayline(irs.left(:,l)',dt,w,conf)';
        brir(:,2) = brir(:,2) + delayline(irs.right(:,l)',dt,w,conf)';
    end

    brir(:,1)=brir(:,1)*10^(Af(ii)/20);
    brir(:,2)=brir(:,2)*10^(-Af(ii)/20);

    % store result to output variable
    brs(:,(ii-1)*2+1:ii*2) = brir;
end

%% ===== Pre-equalization ===============================================
brs = wfs_preequalization(brs,conf);

%% ===== Headphone compensation =========================================
brs = compensate_headphone(brs,conf);
