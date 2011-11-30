% Checks all irs data sets in the current directory, by plotting the relevant
% data

files = dir('*.mat');

for ii = 1:length(files)

    fprintf(1,'\n\nTesting dataset: %s\n',files(ii).name);

    irs = read_irs(files(ii).name);

    % Plot the whole data set
    figure; imagesc(irs.left(1:1500,:)); title('Left Ear');
    figure; imagesc(irs.right(1:1500,:)); title('Right Ear');

    % Plot impulse response
    figure; plot(irs.left(1:1500,100),'b-',irs.right(1:1500,100),'r-');

    % Extract ITD and ILD
    %itd = interaural_time_difference(irs.left,irs.right,irs.fs);
    ild = interaural_level_difference(irs.left,irs.right);
    %plot_itd(itd);
    angle1 = degree(irs.apparent_azimuth(1))
    angle2 = degree(irs.apparent_azimuth(end))
    length(irs.apparent_azimuth)
    %plot_ild(phi,ild);
    figure; plot(degree(irs.apparent_azimuth),ild);

    key_in = input('Resume to next by pressing Enter');
    close all;

end
