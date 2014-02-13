function outsig = auralize_ir(ir,content,usenorm,conf)
%AURALIZE_IR auralizes an impulse response with an audio file/signal
%
%   Usage: outsig = auralize_ir(ir,[[[content],normalize],conf])
%
%   Input parameters:
%       ir        - impulse response (IR). Also an binaural room scanning
%                   (BRS) matrix can be auralized, then only the two first
%                   channels will be used.
%       content   - content file or signal vector to be used for auralisation
%                   (mono, if it contains more than one channel, only the
%                   first will be used).
%                   Also predefined content can be used by applying the
%                   one of the following strings:
%                   'speech', 'noise', 'pinknoise', 'cello', 'castanets'.
%                   Then these contents will be used to auralise the IR.
%                   The corresponding content files are specified in
%                   SFS_config.
%       normalize - normalize the signal (1 or 0), default: 1
%       conf      - optional configuration struct (see SFS_config)
%
%   AURALIZE_IR(ir,content,normalize) convolves the first two channels of the
%   given IR with the given content and returns the resulting outsig. If
%   instead of an explicite content file or vector only a string containig
%   'speech', 'noise', 'pinknoise', 'cello' or 'castanets' is given, the
%   corresponding content file as specified in conf is used.
%
%   see also: auralize_ir_file, ir_wfs, ir_generic, ir_point_source

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
nargmin = 2;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(ir);

if nargin<3
    usenorm = 1;
    conf = SFS_config;
elseif nargin<nargmax
    conf = SFS_config;
end
isargscalar(usenorm);
isargstruct(conf);


%% ===== Configuration ==================================================
% Sampling rate
fs = conf.fs;
% Auralisation files are used directly in the code below in order to made these
% settings not neccessary


%% ===== Get the right content ==========================================
if isnumeric(content)
    contentfs = conf.fs;
else
    if strcmp(content,'castanets')
        contentfile = conf.castanetsfile;
    elseif strcmp(content,'speech')
        contentfile = conf.speechfile;
    elseif strcmp(content,'cello')
        contentfile = conf.cellofile;
    elseif strcmp(content,'noise')
        contentfile = conf.noisefile;
    elseif strcmp(content,'pinknoise')
        contentfile = conf.pinknoisefile;
    elseif ~exist(content,'file')
        error('%s: %s file was not found.',upper(mfilename),content);
    else
        contentfile = content;
    end
    % Read the content file
    [content,contentfs] = wavread(contentfile);
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
