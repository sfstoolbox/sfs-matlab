function x0 = secondary_source_positions(L,conf)
%SECONDARY_SOURCE_POSITIONS Generates the positions and directions of the
%   secondary sources
%   
%   Usage: x0 = secondary_source_positions(L,conf)
%          x0 = secondary_source_positions(L)
%
%   Input options:
%       L           - the size (m) of the array (length for a linear array,
%                     diameter for a circle or a box)
%       conf        - configuration struct
%
%   Output options:
%       x0          - secondary source positions and directions (m)
%
%   SECONDARY_SOURCES_POSITIONS(L) generates the positions and directions x0
%   of secondary sources for a given geometry (conf.array) and array size (L).
%   Alternatively, if conf.x0 is set, it returns the positions and directions
%   specified there.
%   The direction of the sources is given as a point. the direction can be
%   estimated by the vector pointing from the secondary source position point to
%   this direction point.
%
%   Geometry (for the linear array):
%
%                                y-axis
%                                   ^
%                                   |
%                                   |
%                                   |
%          v--v--v--v--v--v--v--v--v|-v--v
%          |              X0        |
%    (Loudspeaker)  (Array center)  |
%                                   |
%       --------------------------------------------------------> x-axis
%
% see also: secondary_source_selection, secondary_source_number, tapwin
%

% AUTHOR: Sascha Spors, Hagen Wierstorf

% NOTE: If you wanted to add a new type of loudspeaker array, do it in a way,
% that the loudspeakers are ordered in a way, that one can go around for closed
% arrays. Otherwise the tapering window function will not work properly.


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(L);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================

% Array type
array = conf.array;
% Center of the array
X0 = position_vector(conf.X0);
% Distance between secondary sources
dx0 = conf.dx0;
% Given secondary sources
x0 = conf.x0;


%% ===== Calculation =====================================================

% Check if we have already predefined secondary sources
if length(x0)>0
    isargsecondarysource(x0);
    % If we have predefined secondary sources leave this function
    return
end

% Get the number of secondary sources
[nls, L] = secondary_source_number(L,conf);

x0 = zeros(nls,6);
if strcmp('linear',array)
    % === Linear array ===
    % Positions of the secondary sources
    x0(:,1) = X0(1) + linspace(-L/2,L/2,nls)';
    x0(:,2) = X0(2) * ones(nls,1);
    x0(:,3) = X0(3) * ones(nls,1);
    % Direction of the secondary sources
    x0(:,4) = x0(:,1);
    x0(:,5) = x0(:,2) + 1;
    x0(:,6) = x0(:,3);
    % === Add jitter to the loudspeaker positions ===
    %x0(1) = X0-size/2;
    %x0(nls) = X0+size/2;
    %for ii=2:nls-1
    %    jitter = size/(4*nls) * randn;
    %    x0(ii) = X0-size/2+(ii-1)*size/nls + jitter;
    %end
elseif strcmp('circle',array)
    % === Circular array ===
    % Positions of the secondary sources
    phi = linspace(0,(2-2/nls)*pi,nls)'; % 0..2pi
    theta = zeros(nls,1);
    [cx,cy,cz] = sph2cart(phi,theta,L/2);
    x0(:,1:3) = [cx,cy,cz] + repmat(X0,nls,1);
    % Direction of the secondary sources
    x0(:,4:6) = repmat(X0,nls,1) .* ones(nls,3);
elseif strcmp('box',array)
    % === Boxed loudspeaker array ===
    % Number of secondary sources per linear array
    % FIXME: can nls/4 be another value than an integer?
    nbox = nls/4;
    % Position and direction of the loudspeakers
    % top
    x0(1:nbox,1) = X0(1) + linspace(-L/2,L/2,nbox)';
    x0(1:nbox,2) = X0(2) + ones(nbox,1) * L/2 + dx0;
    x0(1:nbox,3) = X0(3) + zeros(nbox,1);
    x0(1:nbox,4) = x0(1:nbox,1);
    x0(1:nbox,5) = x0(1:nbox,2) - 1;
    x0(1:nbox,6) = x0(1:nbox,3);
    % right
    x0(nbox+1:2*nbox,1) = X0(1) + ones(nbox,1) * L/2 + dx0;
    x0(nbox+1:2*nbox,2) = X0(2) + linspace(L/2,-L/2,nbox)';
    x0(nbox+1:2*nbox,3) = X0(3) + zeros(nbox,1);
    x0(nbox+1:2*nbox,4) = x0(nbox+1:2*nbox,1) - 1;
    x0(nbox+1:2*nbox,5) = x0(nbox+1:2*nbox,2);
    x0(nbox+1:2*nbox,6) = x0(nbox+1:2*nbox,3);
    % bottom
    x0(2*nbox+1:3*nbox,1) = X0(1) + linspace(L/2,-L/2,nbox)';
    x0(2*nbox+1:3*nbox,2) = X0(2) - ones(nbox,1) * L/2 - dx0;
    x0(2*nbox+1:3*nbox,3) = X0(3) + zeros(nbox,1);
    x0(2*nbox+1:3*nbox,4) = x0(2*nbox+1:3*nbox,1);
    x0(2*nbox+1:3*nbox,5) = x0(2*nbox+1:3*nbox,2) + 1;
    x0(2*nbox+1:3*nbox,6) = x0(2*nbox+1:3*nbox,3);
    % left
    x0(3*nbox+1:nls,1) = X0(1) - ones(nbox,1) * L/2 - dx0;
    x0(3*nbox+1:nls,2) = X0(2) + linspace(-L/2,L/2,nbox)';
    x0(3*nbox+1:nls,3) = X0(3) + zeros(nbox,1);
    x0(3*nbox+1:nls,4) = x0(3*nbox+1:nls,1) + 1;
    x0(3*nbox+1:nls,5) = x0(3*nbox+1:nls,2);
    x0(3*nbox+1:nls,6) = x0(3*nbox+1:nls,3);
elseif strcmp('U',array)
    to_be_implemented(mfilename);
else
    error('%s: %s is not a valid array type.',upper(mfilename),array);
end
