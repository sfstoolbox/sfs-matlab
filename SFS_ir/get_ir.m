function ir = get_ir(irs,phi,delta,r,X0)
%GET_IR returns a IR for the given apparent angle
%
%   Usage: ir = get_ir(irs,phi,[delta,[r]])
%
%   Input parameters:
%       irs     - IR data set
%       phi     - azimuth angle for the desired IR (rad)
%       delta   - elevation angle for the desired IR (rad)
%       r       - distance for the desired IR (m)
%
%   Output parameters:
%       ir      - IR for the given angles (length of IR x 2)
%
%   GET_IR(irs,phi,delta,r) returns a single IR set for the given angles phi and
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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    r = 1;
elseif nargin==nargmax-2
    r = 1;
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
if length(irs.distance) == 1
    irs.distance = repmat(irs.distance,1,length(irs.apparent_azimuth));
end

idx = findrows(...
    roundto([irs.apparent_azimuth' irs.apparent_elevation' irs.distance'],prec),...
    roundto([phi,delta,r],prec));
if idx
    if length(idx)>1
        error(['%s: the irs data set has more than one entry corresponding ',...
               'an azimuth of %.3f deg and an elevation of %.3f deg.'],...
            upper(mfilename),degree(phi),degree(delta));
    end
    ir(:,1) = irs.left(:,idx);
    ir(:,2) = irs.right(:,idx);


else
    
   % get the three nearest IRs via scalar product
   [ir1,ir2,ir3,x0,desired_point] = findnearestneighbour_scalar(irs,phi,delta,r,3,X0);
   L = [x0(:,1),x0(:,2),x0(:,3)];
   
   
   if rank(L) == 3
          
       ir = intpol_ir3d(desired_point,ir1,ir2,ir3,L);
  
   else
          
       ir = intpol_ir2d(ir1,ir2,x0,desired_point);
       
   end

       
end

    % Check if we need to do 2D or 3D interpolation
    % FIXME:
    % Ein problem bleibt hiernoch. Ich weiß nicht genau, hattest Du das nicht
    % schon mal abgefangen? Folgendes Problem:
    % x------------O--x--x--x
    % O ist der gewollte Punkt, x sind die gemessenen Punkte. In diesem Fall
    % würden die ersten drei Nachbarn alle rechts von dem Punkt liegen, was aber
    % Schwachsinn ist für interpolieren, da wir dafür auch den Punkt ganz links
    % benötigen.
%     if all(neighbours(1)-neighbours(2)) && all(neighbours(2)-neighbours(3))
%         % 3D interpolation
%         
%     else
%         to_be_implemented();
%         % FIXME:
%         % An dieser Stelle sollten wir testen, ob die 3 nächsten Nachbarn eine
%         % Ebene aufspannen, die durch den Nullpunkt geht. Wenn dem so ist,
%         % benötigen wir nur eine Interpolation zwischen den dichtesten 2
%         % Punkten.
%         %ir = intpol_ir(desired_point, ...
%         %    x0(:,idx(1)),[irs.left(:,idx(1)) irs.right(:,idx(1))], ...
%         %    x0(:,idx(2)),[irs.left(:,idx(2)) irs.right(:,idx(2))]);
%     end

end % of main function

%% ===== Subfunctions ====================================================
% round the input matrix m to the given precission prec in degree
function m = roundto(m,prec)
    m = round(degree(m)/prec)*prec;
end
