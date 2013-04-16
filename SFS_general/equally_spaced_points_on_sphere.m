function x0 = equally_spaced_points_on_sphere(L,conf)

% x0 = equally_spaced_points_on_sphere(conf) loads the coordinates and
% weights of N equally spaced points on a sphere. N has to be set in conf
% (conf.number_of_points_on_sphere).

% Inputs: L: diameter of spherical array
%         conf.number_of_points_on_sphere: has to be a sqaure number in the range
%                                          of 1...6561
%         conf.X0: reference point

% Output: (Nx7)-Matrix containing x-,y- and z-coordinates of N equally
%          spaced points on a sphere with diameter L in its first 3 columns and
%          its weights in the 7th column. The columns 4-6 hold the x-,y-,z-coordinates
%          of the appropriate normal vector with respect to the columns
%          1-3.

%% error check
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));

if mod(conf.number_of_points_on_sphere,sqrt(conf.number_of_points_on_sphere)) ~= 0
    error('%s: conf.number_of_points_on_sphere has to be a squared number.',upper(mfilename));

    elseif isempty(L)    
    error('%s: diameter of the desired sphere has to be set.',upper(mfilename));
    
else
%% re-define conf.variables
N = conf.number_of_points_on_sphere;
%% load file and calculate x0
% get current path
basepath = which('equally_spaced_points_on_sphere' );
basepath = basepath(1:end-33);
addpath([basepath,'MinimumEnergyPointsOnASphere']);
N = [num2str(N) 'points.mat'];
x0 = load(N,'-ascii');
x0 = L/2*x0;
% x0(:,1:3) = L/2*x0(:,1:3);
x0(:,7) = x0(:,4);
x0(:,8) = sin(acos(x0(:,3)./(L/2)));
%% get normal vectors of x0
X0 = conf.X0;
x0(:,4:6) = direction_vector(x0(:,1:3),repmat(X0,length(x0),1).*ones(length(x0),3)); 
end

end % end of function
