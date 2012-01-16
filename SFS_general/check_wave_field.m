function check_wave_field(P,frame)
%CHECK_WAVE_FIELD checks if we have any activity in the wave field and returns a
%   warning otherwise.
%
%   Usage: check_wave_field(P,frame)
%
%   Input parameters:
%       P       - wave field
%       frame   - used time frame
%
%
%   CHECK_WAVE_FIELD(P,frame) checks if the wave field is different from zero.
%   If this is not the case it returns a warning.
%
%   see also: wave_field_imp_wfs_25d, norm_wave_field
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(P);
isargscalar(frame);


%% ===== Computation =====================================================
if max(abs(P(:)))==0 | isnan(P(:))
    warning('SFS:check_wave_field',...
        ['The activity in the simulated wave field is zero. ',...
         'Maybe you should use another time frame than %i. ', ...
         'You can set the time frame with conf.frame.'],frame);
end
