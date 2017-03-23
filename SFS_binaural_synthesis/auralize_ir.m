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
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
    [content,contentfs] = audioread(contentfile);
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
