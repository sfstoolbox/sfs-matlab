function check_wave_field(P)
%CHECK_WAVE_FIELD checks if we have any activity in the wave field
%   Usage: bool = check_wave_field(P,x,y,yref)
%
%   Input parameters:
%       P       - wave field
%
%
%   CHECK_WAVE_FIELD(P) checks if the wave field is different from zero. If this
%   is not the case it returns a warning.
%
%   see also: wave_field_imp_wfs_25d, norm_wave_field

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(P);


%% ===== Computation =====================================================
if max(abs(P(:)))==0
    warning('SFS:check_wave_field',...
        ['The activity in the simulated wave field is zero. ',...
         'Maybe you should use another time frame conf.frame.']);
end
