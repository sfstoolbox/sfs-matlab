function ir = get_ir(irs,phi,delta)
%GET_IR returns a IR for the given apparent angle
%
%   Usage: ir = get_ir(irs,phi,[delta])
%
%   Input parameters:
%       irs     - IR data set
%       phi     - azimuth angle for the desired IR (rad)
%       delta   - elevation angle for the desired IR (rad)
%
%   Output parameters:
%       ir      - IR for the given angles (length of IR x 2)
%
%   GET_IR(irs,phi,delta) returns a single IR set for the given angles phi and
%   delta. If the desired angles are not present in the IR data set an
%   interpolation is applied to create the desired angles.
%
%   see also: read_irs, slice_irs, ir_intpol

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
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin))
if nargin==nargmax-1
    delta = 0;
end


%% ===== Configuration ==================================================
% Precission of the wanted angle. If a IR within the given precission could be
% found no interpolation is applied.
prec = 0.1; % degree


%% ===== Computation ====================================================

% === Check the given angles ===
% Ensure -pi <= phi < pi and -pi/2 <= delta <= pi/2
phi = correct_azimuth(phi);
delta = correct_elevation(delta);

% === IR interpolation ===
% Check if the IR dataset contains a measurement for the given angles
% phi and delta. If this is not the case, interpolate the dataset for the given
% angles.

% If azimuth and elevation could be found
idx = findrows(...
    roundto([irs.apparent_azimuth' irs.apparent_elevation'],prec),...
    roundto([phi,delta],prec));
if idx
    if length(idx)>1
        error(['%s: the irs data set has more than one entry corresponding ',...
               'an azimuth of %.3f deg and an elevation of %.3f deg.'],...
            upper(mfilename),degree(phi),degree(delta));
    end
    ir(:,1) = irs.left(:,idx);
    ir(:,2) = irs.right(:,idx);


else (isempty(idx));
    
    % calculate the x-,y-,z-coordinates for the desired angles(phi,delta)
    [x,y,z] = sph2cart(phi,delta,1);
    vec = [x,y,z];
    
    % calculate the x-,y-,z-coordinates for all known IRs
    [x,y,z] = sph2cart(irs.apparent_azimuth(1,:),irs.apparent_elevation(1,:),irs.distance);
    x0 = [x;y;z];
    
    % get the maxima in order to get the indices of the three closest IRs 
    % with respect to the desired angles(phi,delta)
    scalar_product = vec*x0;
    [~,index1] = max(scalar_product);
    if length(index1)>1
        index1 = index1(1,1);
    end
    
    scalar_product(:,index1) = -1-eps;
    [~,index2] = max(scalar_product);
    if length(index2)>1
        index2 = index2(1,1);
    end
       
    scalar_product(:,index2) = -1-eps;
    [~,index3] = max(scalar_product);
    if length(index3)>1
        index3 = index3(1,1);
    end
       
fprintf('The indices of the found IRs with get_ir_new_version are: index1 = %d, index2 = %d,index3 = %d\n',index1,index2,index3);     
        
    % get the desired IRs
    ir1(:,1) = irs.left(:,index1);
    ir1(:,2) = irs.right(:,index1);

    ir2(:,1) = irs.left(:,index2);
    ir2(:,2) = irs.right(:,index2);

    ir3(:,1) = irs.left(:,index3);
    ir3(:,2) = irs.right(:,index3);
           
    % IR interpolation
    ir = intpol_ir(ir1,irs.apparent_azimuth(index1),irs.apparent_elevation(index1),...
                   ir2,irs.apparent_azimuth(index2),irs.apparent_elevation(index2),...
                   ir3,irs.apparent_azimuth(index3),irs.apparent_elevation(index3),...    
                   phi,delta...
                   );
end
end % of main function

%% ===== Subfunctions ====================================================
% round the input matrix m to the given precission prec in degree
function m = roundto(m,prec)
    m = round(degree(m)/prec)*prec;
end