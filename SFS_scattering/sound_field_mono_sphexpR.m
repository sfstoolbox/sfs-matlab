function [P, x, y, z] = sound_field_mono_sphexpR(X,Y,Z,Al,f,x0,conf)


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z,Al);
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
[Jn, ~, Ynm, x, y, z] = sphbasisXYZ(X,Y,Z,f,x0,conf);

P = sound_field_mono_sphbasis(Al,Jn,Ynm,conf);
end

