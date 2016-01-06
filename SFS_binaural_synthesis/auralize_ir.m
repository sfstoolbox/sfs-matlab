function outsig = auralize_ir(ir,content,usenorm,conf)
%AURALIZE_IR auralizes an impulse response with an audio file/signal
%
%   Usage: outsig = auralize_ir(ir,content,[normalize],conf)
%
%   Input parameters:
%       ir        - impulse response (IR). Also an binaural room scanning
%                   (BRS) matrix can be auralized, then only the two first
%                   channels will be used.
%       content   - content file or signal vector to be used for auralisation
%                   (mono, if it contains more than one channel, only the
%                   first will be used).
%       normalize - normalize the signal (1 or 0), default: 1
%       conf      - configuration struct (see SFS_config)
%
%   AURALIZE_IR(ir,content,normalize,conf) convolves the first two channels of
%   the given IR with the given content and returns the resulting outsig. The
%   content can be specified in the form of an explicite content file or vector.
%
%   See also: ir_wfs, ir_generic, ir_point_source

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input parameters and configuration =================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(ir);
if nargin<nargmax
    conf = usenorm;
    usenorm = 1;
end
isargscalar(usenorm);
isargstruct(conf);


%% ===== Configuration ==================================================
% Sampling rate
fs = conf.fs;


%% ===== Get the right content ==========================================
if isnumeric(content)
    contentfs = conf.fs;
elseif ~exist(content,'file')
    [content,contentfs] = wavread(contentfile);
else
    error('%s: %s file was not found.',upper(mfilename),content);
end


%% ===== Convolution of the IR with content =============================
% Check if we have to resample
if contentfs~=fs
    content = resample(content,fs,contentfs);
end
% Check if the content is only one channel and fix it otherwise
if min(size(content))~=1
    if size(content,1)<size(content,2)
        content = content(1,:);
    else
        content = content(:,1);
    end
end
% Convolve the two
if size(ir,2)>2
    warning('Your impulse response has more than two channels.');
end
% Convolve the impulse responses with the content signal
outsig = convolution(ir,content);
% Scale output
if(usenorm)
    outsig = norm_signal(outsig);
end
