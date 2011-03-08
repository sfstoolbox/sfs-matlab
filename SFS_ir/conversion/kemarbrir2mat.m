function kemarbrir2mat(irsset,irspath)
%KEMARBRIR2MAT converts IRs data given by our KEMAR to our mat-file based format
%
%   Usage: kemarbrir2mat(irsset,irspath);
%
%   Input options:
%       irsset  - IR sets measured with VariSphear and KEMAR. Currently the 
%                 following are available:
%                   'spirit1_src1'  - BRIR of Spirit with source at [0 2]
%                   'spirit1_src2'  - BRIR of Spirit with source at [-sqrt(3) 2]
%                   'spirit1_src3'  - BRIR of Spirit with source at [sqrt(3) 2]
%                 NOTE: you still have to give the matching path to the given
%                 data set!
%       irspath - path to the directory containing the IR data
%
%   KEMARBRIR2MAT(irsset,irspath) converts the IRs data given by the irsset
%   and stored at the given irspath in our own mat-file based format. See:
%   https://dev.qu.tu-berlin.de/projects/sfs/wiki/IRs_mat-file_format

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargchar(irsset);
isargdir(irspath);


%% ===== Computation =====================================================
% Dir to save the IRs
outdir = 'ir_databases';

if strcmp(irsset,'spirit1_src1')
    irs.description = ...
        ['KEMAR measurement in Spirit of T-Labs Berlin. Used elevation: ', ...
         '0; azimuth: -90:1:90. Rotation: head. Ears: large. Source at ',...
         '0 degree and 2 m'];
    irfilebase = '2011-02-07_15-27-40_src1';
elseif strcmp(irsset,'spirit1_src2')
    irs.description = ...
        ['KEMAR measurement in Spirit of T-Labs Berlin. Used elevation: ', ...
         '0; azimuth: -90:1:90. Rotation: head. Ears: large. Source at ',...
         '30 degree and 2 m'];
    irfilebase = '2011-02-07_15-27-40_src2';
elseif strcmp(irsset,'spirit1_src3')
    irs.description = ...
        ['KEMAR measurement in Spirit of T-Labs Berlin. Used elevation: ', ...
         '0; azimuth: -90:1:90. Rotation: head. Ears: large. Source at ',...
         '-30 degree and 2 m'];
    irfilebase = '2011-02-07_15-27-40_src3';
else
    error('%s: the given irsset is not available.',upper(mfilename));
end

% Read one angle to get the independent variables
irfile = sprintf('%s/%s_head+000.mat',irspath,irfilebase);
load(irfile);

irs = data;
irs = rmfield(irs,'ir');
irs = rmfield(irs,'ir_info');
% Read the data
for ii = 1:181

    if ii<91
        irfile = sprintf('%s/%s_head-%03.0f.mat',irspath,irfilebase,abs(ii-91));
    else
        irfile = sprintf('%s/%s_head+%03.0f.mat',irspath,irfilebase,ii-91);
    end
    load(irfile);
    irs.apparent_azimuth(ii) = correct_azimuth(data.apparent_azimuth);
    irs.apparent_elevation(ii) = correct_elevation(data.apparent_elevation);
    irs.head_azimuth(ii) = data.head_azimuth;
    irs.left(:,ii) = data.ir(:,1);
    irs.right(:,ii) = data.ir(:,2);

end

% Reorder fields
irs = order_irs_fields(irs);
% Reorder entries
irs = correct_irs_angle_order(irs);
% Check irs format
check_irs(irs);

% Create the outdir
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Write IR mat-file
outfile = sprintf('%s/KEMAR_%s.mat',outdir,irsset);
save('-v7',outfile,'irs');
