function savewav(y,filename,fs)
%SAVEWAV saves floating point audio data in a *.wav with 32bit precision
%
%   Usage: savewav(y,filename,fs)
%
%   Input parameters:
%       y         - audio data in floating point [NFrames x NChannels] 
%       filename  - name of *.wav file (file extension will NOT be added
%                   automatically)
%       fs        - sample rate / Hz 
%
%   Since wavwrite has been removed in MATLAB 2015b and audiowrite does not
%   support audio data with a channel number higher than 256, this is an 
%   alternative to save massive multichannel. This function does only support
%   floating point data and saves it 32bit floating point precision.
%   
%   See also: audiowrite, audioread, fwrite

% Copyright (c) 2015 Fiete Winter
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%% Input Check
narginchk(3,3);
if ~isfloat(y)
  error('%s: only floating point supported', upper(mfilename));
end
if ~ndims(y) > 2
  error('%s: only audio data with [NFrames x NChannels] supported', ...
    upper(mfilename));
end

%%
[NFrames, NChannels] = size(y);  % Number of Frames and Channels
NSamples = NFrames*NChannels;  % total Samples
Nbits = 32;  % 32 floating point
Nbyte = 4;  % bytes per Samples

fileID = fopen(filename,'w', 'ieee-le');  % open wav-file (little Endian)

fwrite(fileID, 'RIFF', 	'char*1');  % RIFF header
% number of bytes following (chunk size)
Noverallbytes = 4 + 4; % for 'WAVE' and 'fmt '
Noverallbytes = Noverallbytes + 4;  % for the number of bytes in 'fmt '
Noverallbytes = Noverallbytes + 16; % bytes in 'fmt '
Noverallbytes = Noverallbytes + 12; % for 'fact'
Noverallbytes = Noverallbytes + 8 + Nbyte*NSamples;  % 'data' chunk
fwrite(fileID, Noverallbytes, 'uint32');  

fwrite(fileID, 'WAVE', 	'char*1');  % WAVE header

%% fmt chunk
fwrite(fileID, 'fmt ', 	'char*1');  % mind the whitespace!
fwrite(fileID, 16, 	'uint32');  % number of bytes following (chunk size)
fwrite(fileID, 3, 'uint16');  % data type (3 for float)
fwrite(fileID, NChannels, 'uint16');  % number of channels
fwrite(fileID, fs, 'uint32');  % samplerate
fwrite(fileID, fs * NChannels * Nbyte, 'uint32');  % bytes per second
fwrite(fileID, NChannels * Nbyte, 'uint16');  % frame size
fwrite(fileID, Nbits, 'uint16');  % bits per sample

%% fact chunk (floating point only)
fwrite(fileID, 'fact', 	'char*1');
fwrite(fileID, 4, 	'uint32');  % number of bytes following (chunk size)
fwrite(fileID, NFrames, 	'uint32');  % number of frames

%% data chunk
fwrite(fileID, 'data', 	'char*1');
fwrite(fileID, Nbyte*NSamples, 'uint32');  % number of bytes following (chunk size)
y = y.';  % prepare data for interlacing
fwrite(fileID, y(:), 'float32');  % interlaced data

%%
fclose(fileID);

end