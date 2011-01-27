function plot_itd(itd,phi)
%PLOT_ITD plots the given ITD values
%   Usage: plot_itd(itd,phi)
%          plot_itd(itd)
%
%   Input options:
%       itd -   vector with given ITD values (e.g. crfeated with extract_itd)
%       phi -   corresponding angles (optional, default: -180°..179°)
%
%   PLOT_ITD(itd,phi) creates a figure with the given ITD values and add
%   corresponding labels etc.
%
%   see also: extract_itd, plot_ild

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if ~isnumeric(itd) || ~isvector(itd)
    error('%s: itd has to be a vector!',upper(mfilename));
end
if nargin==nargmax
    if ~isnumeric(phi) || ~isvector(phi)
        error('%s: phi has to be a vector!',upper(mfilename));
    end
    if length(itd)~=length(phi)
        error('%s: phi has to have the same length as itd!',upper(mfilename));
    end
else
    phi = -180:1:179;
end


%% ===== Plotting ========================================================
figure;
plot(phi,1000*itd);
axis([phi(1),phi(end),-1,1]);
xlabel('phi (°)');
ylabel('ITD (ms)')
