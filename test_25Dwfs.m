clear variables
close all

conf = SFS_config;
conf.plot.usedb = true;
conf.plot.useplot = false;
conf.usenormalisation = false;

conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.size = 100;
conf.secondary_sources.number = 1000;
conf.secondary_sources.center = [0,5,0];
conf.xref = [0,0,0];

conf.usetapwin = true;
conf.tapwinlen = 0.2;

%%
f = 1000;  % temporal frequency
positions = { [0,7,0], [-sqrt(0.5), -sqrt(0.5), 0], [0,3,0,0,-1,0] };  % source positions
sources = {'ps', 'pw', 'fs'};
gtsources = {'ps', 'pw', 'ps'};

X = [-2,2];
Y = [-2,2];
Z = 0;

for idx=1:length(positions)
  xs = positions{idx};
  src = sources{idx};
  gt = gtsources{idx};
  
  for driving_functions = {'reference_point', 'reference_line'}
    conf.driving_functions = driving_functions{:};

    [Pgt,x,y] = sound_field_mono(X,Y,Z,[xs(1:3),0,-1,0,1], gt, 1, f, conf);
    Pwfs = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
    
    figure;
    imagesc(Y,X,db(1 - Pwfs./Pgt));    
    title(sprintf('%s %s', src, driving_functions{:}), 'Interpreter', 'none');
    set(gca, 'YDir', 'normal');
    
    % plot_sound_field(Pwfs./Pgt, X, Y, Z, [], conf);   
  end 
end