function [J, H2, Y, x, y, z]  = eval_sphbasis_mono_XYZgrid(X,Y,Z,f,xq,conf)


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  xq = [0, 0, 0];
end  
isargposition(xq);

%% ===== Computation ====================================================

% Create a x-y-grid
[xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,conf);

k = 2*pi*f/conf.c;  % wavenumber

% shift coordinates to expansion coordinate
xx = xx-xq(1);
yy = yy-xq(2);
zz = zz-xq(3);

% coordinate transformation
r = vector_norm(cat(3,xx,yy,zz),3);
phi = atan2(yy,xx);
theta = asin(zz./r);

[J, H2, Y] = eval_sphbasis_mono(r, theta, phi, k, conf);

end




