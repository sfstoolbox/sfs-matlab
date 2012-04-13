function ir = fix_ir_length(ir,N,dt)
%FIX_IR_LENGTH pads zeros or removes entries from the IR according to length N
%
%   Usage: ir = fix_ir_length(ir,N,dt)
%
%   Input parameters:
%       ir  - impulse response (IR)
%       N   - number of samples the calculated BRIR should have
%       dt  - time delay for the given setup the IR will be shifted with
%
%   Output paramteres:
%       ir  - corrected IR
%
%   FIX_IR_LENGTH(IR,N,DT) pads zeros or removes the end of the given IR in
%   order to have a IR with the correct length to feed it with the time delay dt
%   in the desired BRIR with the length N.
%
%   see also: brs_point_source, get_ir

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Fix IR ==========================================================
% length of IR
lenir = length(ir(:,1));
% append zeros if to short
if(lenir<N-dt)
    ir = cat(1,ir,zeros(N-lenir,2));
% remove the end of the IR, if to long
else
    ir=ir(1:N,:);
end
