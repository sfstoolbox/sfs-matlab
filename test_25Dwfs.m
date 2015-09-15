clear variables

conf = SFS_config_example;
conf.plot.usedb = true;
conf.plot.useplot = false;
conf.usenormalisation = false;

conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.size = 100;
conf.secondary_sources.number = 5000;
conf.secondary_sources.center = [0,1,0];

conf.usetapwin = false;

%%
f = 5000;  % temporal frequency
positions = { [0,3,0], [0, -1, 0] };  % source positions
sources = {'ps', 'pw'};

for idx=1:2
  xs = positions{idx};
  src = sources{idx};
  
  x0 = secondary_source_positions(conf);
  x0 = secondary_source_selection(x0, xs, src);
  x0 = secondary_source_tapering(x0, conf);  
  
  for driving_functions = {'reference_point', 'reference_line'}
    conf.driving_functions = driving_functions{:};

    D0 = driving_function_mono_wfs(x0,xs,src,f,conf);  

    X = [-2,2];
    Y = 0;
    Z = 0;

    [Pgt,x] = sound_field_mono(X,Y,Z,[xs,0,-1,0,1], src, 1, f, conf);
    Pwfs = sound_field_mono(X,Y,Z,x0, 'ps', D0, f, conf);

    figure;
    subplot(1,2,1)
    plot(x, db(Pwfs), 'r');
    hold on;
    plot(x, db(Pgt), 'b');
    hold off;

    X = 0;
    Y = [-1,3];
    Z = 0;  

    [Pgt,~,y] = sound_field_mono(X,Y,Z,[xs,0,-1,0,1], src, 1, f, conf);
    Pwfs = sound_field_mono(X,Y,Z,x0, 'ps', D0, f, conf);

    subplot(1,2,2)
    plot(y, db(Pwfs), 'r');
    hold on;
    plot(y, db(Pgt), 'b');
    hold off;
    
  end
  
end