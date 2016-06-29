clear variables
close all

conf = SFS_config;
conf.plot.usedb = true;
conf.plot.useplot = false;
conf.usenormalisation = false;

conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.size = 20;
conf.secondary_sources.number = 512;
conf.secondary_sources.center = [0,6,0];
conf.xref = [0,0,0];

conf.usetapwin = true;
conf.tapwinlen = 0.2;

%%
f = 2000;  % temporal frequency
positions = { [0,9,0], [0,-1,0], [0,3,0,0,-1,0] };  % source positions
sources = {'ps', 'pw', 'fs'};
gtsources = {'ps', 'pw', 'ps'};

X = [-2,2];
Y = [-2,2];
Z = 0;

for idx=1:length(positions)
  xs = positions{idx};
  src = sources{idx};
  gt = gtsources{idx};
  
  figure;
  ddx = 1;
  for driving_functions = {'reference_point', 'reference_line'}
    conf.driving_functions = driving_functions{:};

    [Pgt,x,y] = sound_field_mono(X,Y,Z,[xs(1:3),0,-1,0,1], gt, 1, f, conf);
    Pwfs = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
    
    subplot(1,2,ddx);
    imagesc(Y,X,db(1 - Pwfs./Pgt));    
    title(sprintf('%s %s', src, driving_functions{:}), 'Interpreter', 'none');
    set(gca, 'YDir', 'normal');
    colorbar;
    hold on;
    if strcmp( conf.driving_functions, 'reference_point' )
       plot(conf.xref(1),conf.xref(2),'gx');
    else
       plot(conf.xref(1)+X,conf.xref([2,2]),'g--');
    end    
    hold off;
    
    ddx= ddx+1;
  end     
end