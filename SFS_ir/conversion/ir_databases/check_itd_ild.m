% Check the ITD and ILD of the data bases in order to check the orientation of
% the KEMAR manikin
itd0 = [];
ild0 = [];
files = dir('*.mat');

% Extract ITD and ILD
for ii = 1:length(files)
    irs = read_irs(files(ii).name);
    itd = interaural_time_difference(irs.left,irs.right,irs.fs,'hilbert');
    ild = interaural_level_difference(irs.left,irs.right);
    idx = find(itd<0,30,'first');
    for jj = 1:30
        % Check if we have a negative slope to ensure we are really near a value
        % of 0 degree
        if itd(idx(jj))>itd(idx(jj)+1)
            itd0(ii) = degree(irs.apparent_azimuth(idx(jj)));
            break;
        end
    end
    %itd0(ii) = degree(irs.apparent_azimuth(find(itd(idx-20:idx+20)<0,1,'first')));
    ild0(ii) = ild(find(irs.apparent_azimuth>0,1,'first'));
    fprintf(1,'\n%s \n angle offset: %.2f degree \t ILD: %.2f dB\n', ...
        files(ii).name,itd0(ii),ild0(ii))
end

itd0
ild0
fprintf(1,'\nMean values \n angle offset: %.2f degree \t ILD: %.2f dB\n', ...
    mean(itd0),mean(ild0))
