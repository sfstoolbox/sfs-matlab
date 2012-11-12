conf = SFS_config;
L = 3;
conf.array = 'circle';
X = [0 0];
phi = pi/2;
src = 'ps';
xs = [0 2.5];
conf.xref = [0 0];
conf.usehpre = 1;
conf.usehcomp = 0;
% get dirac impulses as impulse responses
irs = dummy_irs;
% 1. loudspeaker spacing = 0.16m
conf.dx0 = 0.16;
conf.hprefhigh = aliasing_frequency(conf.dx0);
ir = ir_wfs_25d(X,phi,xs,src,L,irs,conf);
[a,p,f] = easyfft(ir(:,1),conf);
figure; semilogx(f,db(a));
% 2. loudspeaker spacing = 0.32m
conf.dx0 = 0.32;
conf.hprefhigh = aliasing_frequency(conf.dx0);
ir = ir_wfs_25d(X,phi,xs,src,L,irs,conf);
[a,p,f] = easyfft(ir(:,1),conf);
figure; semilogx(f,db(a));
% if you want to listen to the result, you should include the headphone
% compensation
conf.usehcomp = 1;
ir = ir_wfs_25d(X,phi,xs,src,L,irs,conf);
% make sure you have set conf.speechfile in SFS_config.m
sig = auralize_ir(ir,'speech');
sound(sig,conf.fs);
