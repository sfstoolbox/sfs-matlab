function [P, x, y, z] = sound_field_mono_sht(X,Y,Z,Dnm,f,conf)
%SOUND_FIELD_MONO_SHT simulates a sound field given with the spherical 
%harmonics transform of the driving function
%
%   Usage: [P, x, y, z] = sound_field_mono_sht(X,Y,Z,Dnm,f,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       Dnm         - spherical harmonics transform of nfchoa driving function
%       f           - frequency in Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - resulting soundfield
%
%   SOUND_FIELD_MONO_SHT(X,Y,Z,Dnm,f,conf)
%
%   see also: sphbasis_mono_grid, sound_field_mono_sphexp

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(Dnm);
isargsquaredinteger(length(Dnm));
isargnumeric(X,Y,Z);
isargpositivescalar(f);
isargstruct(conf);

%% ===== Configuration ==================================================
Xc = conf.secondary_sources.center;
r0 = conf.secondary_sources.size / 2;
useplot = conf.plot.useplot;

%% ===== Computation ====================================================
[x,y,z] = xyz_grid(X,Y,Z,conf);
% find coordinates, which are inside and outside the loudspeaker array
select = sqrt((x-Xc(1)).^2 + (y-Xc(2)).^2 + (z-Xc(3)).^2) <= r0;

if (numel(x) == 1) x = repmat(x, size(select)); end
if (numel(y) == 1) y = repmat(y, size(select)); end
if (numel(z) == 1) z = repmat(z, size(select)); end

P = zeros(size(x));

% regular (interior) domain inside the loudspeaker array
if any(select(:))
  Pnm = sphexp_mono_sht(Dnm,'R',f,conf);
  P(select) = sound_field_mono_sphexp(x(select), y(select), z(select), ...
    Pnm, 'R', f, Xc,conf);
end

% singular (exterior) domain outside the loudspeaker array
if any(~select(:))
  Pnm = sphexp_mono_sht(Dnm,'S',f,conf);
  if sum( ~select(:) ) == 2
    % this handle cases, where x(~select) would only contain 2 entries, which
    % would be interpreted as a range by the sound_field_... function
    xtmp = x(~select);
    ytmp = y(~select);
    ztmp = z(~select);
    Ptmp(1) = sound_field_mono_sphexp(xtmp(1), ytmp(1), ztmp(1), ...
      Pnm, 'S', f, Xc,conf);
    Ptmp(2) = sound_field_mono_sphexp(xtmp(2), ytmp(2), ztmp(2), ...
      Pnm, 'S', f, Xc,conf);
    P(~select) = [Ptmp(1); Ptmp(2)];
  else
    P(~select) = sound_field_mono_sphexp(x(~select), y(~select), z(~select), ...
      Pnm, 'S', f, Xc,conf);
  end
end

%% ===== Plotting =======================================================
if (nargout==0 || useplot)
    plot_sound_field(P,X,Y,Z,[],conf);
end

end