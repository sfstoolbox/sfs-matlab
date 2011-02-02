function irs = dummy_irs()
% DUMMY_IRS creates a dummy dirac pulse IR set
%   Usage: irs = dummy_irs()
%
%   Output parameters:
%       irs   - irs struct
%
%   DUMMY_IRS() creates a dummy IR data set (Dirac impulse) to check
%   processing without IRs.
%
%   See also: new_irs, IR_format.txt

%   AUTHOR: Hagen Wierstorf


%% ===== Computation =====================================================

% Create dirac pulse
nsamples = 1024;
ir = zeros(nsamples,1);
ir(300) = 1;

irs = new_irs();
for ii=0:360
    irs.left(:,ii+1) = ir;
    irs.right(:,ii+1) = ir;
    irs.apparent_azimuth(ii+1) = correct_azimuth(ii/180*pi);
    irs.apparent_elevation(ii+1) = correct_elevation(0);
end
irs.description = ['HRIR dummy set (Dirac pulse) for testing your',...
                   'frequency response, etc.'];
% Reorder entries
irs = correct_irs_angle_order(irs);
