function [phi,delta] = correct_angle(phi)
%CORRECT_ANGLE is an obsolete function
%
%   It was replaced by correct_azimuth and correct_elevation. use these
%   functions instead!

% AUTHOR: Hagen Wierstorf

error(['%s: Obsolete function! Use correct_azimuth or correct_elevation ',...
       'instead!'],upper(mfilename));
