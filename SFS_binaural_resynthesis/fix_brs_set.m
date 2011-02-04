function fix_brs_set(brsfile) 
%FIX_BRS Fixes a BRS set with a wrong angle direction
%   Usage: fix_brs(brsfile)
%
%   Input parameters:
%       brsfile - BRS file name to fix
%
%   FIX_BRS_SET(brsfile) fixes a BRS file with a wrong angle (clockwise) 
%   direction. This means if a source is moving when you uses the headtracker 
%   but the source shouldn't move. Hopefully this script fixes this 
%   behavior.

% see also: ref_brs_set, wfs_brs_set
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargfile({brsfile},{'brsfile'});


%% ===== Computation =====================================================

% Read brsfile size
[nchannels, nsamples] = wavread (brsfile, 'size');
% Check if we have the right number of channels
if nchannels~=720
    error('%s: The number of channels in the brsfile is %i and not 720!', ...
        upper(mfilename),nchannels);
end
% Read brsfile
[brs, fs, bits] = wavread (brsfile);

% Switch the channels
brs_new = zeros(size(brs));
brs_new(:,1:2) = brs(:,1:2);
for i = 3:2:719
    brs_new(:,i:i+1) = brs(:,722-i:723-i);
end
% NOTE: no new scaling is requiered, because the first two channels are the 
% same as in the old BRS file!

% Write new BRS file
outfile = ['fixed_',brsfile];
wavwrite (brs_new,fs,bits,outfile)

