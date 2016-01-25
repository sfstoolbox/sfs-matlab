function boolean = test_wfs_iir_prefilter()
%TEST_WFS_IIR_PREFILTER tests the IIR WFS pre-equalization filter
%
%   Usage: boolean = test_wfs_iir_prefilter()
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_WFS_IIR_PREFILTER() test the WFS pre-euqalization IIR filter
%   design. This works only in Matlab as the Signal Processing Toolbox is used.
%   See wfs_iir_prefilter.m for details

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



%% ===== Configuration ===================================================
boolean = false;
conf = SFS_config;


%% ===== Calculation =====================================================
% call with default values
hpre1 = wfs_iir_prefilter(conf)
conf.fs = 44100;
conf.hpreflow = 200;
conf.hprefhigh = 1500;
conf.hpreBandwidth_in_Oct = 2;
conf.hpreIIRorder = 4;
hpre2 = wfs_iir_prefilter(conf)
boolean = true;
