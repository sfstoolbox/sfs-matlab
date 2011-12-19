% Check the ITD and ILD of the data bases in order to check the orientation of
% the KEMAR manikin
itd0 = [];
ild0 = [];
files = dir('*.mat');

% Extract ILDs
for ii = 1:length(files)
    irs = read_irs(files(ii).name);
    ild = interaural_level_difference(irs.left,irs.right);
    ild0(ii) = ild(find(irs.apparent_azimuth>0,1,'first'));
end

% Correct ILD and save irs sets
ild_dev = mean(ild0)
%ild_dev = -2.00;
for ii = 1:length(files)
    irs = read_irs(files(ii).name);
    irs.right = gaindb(irs.right,-ild_dev);
    save_irs(irs,files(ii).name);
end
