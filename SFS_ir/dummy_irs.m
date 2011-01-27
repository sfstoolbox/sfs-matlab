function irs = dummy_irs()
% DUMMY_IR Creates a dummy IR set
%   Usage: irs = dummy_irs()
%
%   Output parameters:
%       irs   - 1 x 360 struct containing .left, .right, .angle, .r0, .tag,
%               .description
%
%   DUMMY_IRS() creates a dummy IR data set (Dirac impulse) to check
%   processing without IRs.
%
%   See also: read_irs, create_irs_mat
%

%   AUTHOR: Hagen Wierstorf


%% ===== Computation =====================================================

nsamples = 11025;
ir = zeros(nsamples,1);
ir(300) = 1;

irs = {};
for a=0:359
    irs.left(:,a+1) = ir;
    irs.right(:,a+1) = ir;
    irs.angle(:,a+1) = [a;0];
end
% Has this value any effect
irs.r0 = 1;
irs.tag = 'HRIR';
irs.description = ['HRIR dummy set (Dirac pulse) for testing your',...
                   'frequency response, etc.'];
