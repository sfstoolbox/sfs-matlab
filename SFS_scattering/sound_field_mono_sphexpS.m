function [P, x, y, z] = sound_field_mono_sphexpS(X,Y,Z,Bl,f,x0,conf)


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z,Bl);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  x0 = [0, 0, 0];
end  
isargposition(x0);

%% ===== Computation ====================================================
[~, Hn, Ynm, x, y, z] = sphbasis_mono_XYZgrid(X,Y,Z,f,x0,conf);

P = sound_field_mono_sphbasis(Bl,Hn,Ynm,conf);
end

