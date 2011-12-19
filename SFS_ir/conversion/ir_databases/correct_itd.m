% Check the ITD and ILD of the data bases in order to check the orientation of
% the KEMAR manikin
itd0 = [];
ild0 = [];
files = dir('*.mat');


% Extract ITD
for ii = 1:length(files)
    irs = read_irs(files(ii).name);
    itd = interaural_time_difference(irs.left,irs.right,irs.fs,'hilbert');
    idx = find(itd<0,30,'first');
    for jj = 1:30
        % Check if we have a negative slope to ensure we are really near a value
        % of 0 degree
        if itd(idx(jj))>itd(idx(jj)+1)
            itd0(ii) = degree(irs.apparent_azimuth(idx(jj)));
            break;
        end
    end
end

% Correct ITD
itd_dev = rad(mean(itd0));
for ii = 1:length(files)
    irs = read_irs(files(ii).name);
    irs.apparent_azimuth = correct_azimuth(irs.apparent_azimuth-itd_dev);
    irs = correct_irs_angle_order(irs);
    save_irs(irs,files(ii).name);
end
