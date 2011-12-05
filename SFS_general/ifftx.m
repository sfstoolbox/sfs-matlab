function [ out ] = ifftx( varargin )
%IFFTX Inverse spatial fourier transform
%   Detailed explanation goes here

% FIXME: maybe this function is obsolete?

if nargin == 1
    input_data = varargin{1};
    points     = [];
    dimension  = 1;

elseif nargin == 2
    input_data = varargin{1};
    points     = varargin{2};
    dimension  = 1;
    
elseif nargin == 3
    input_data = varargin{1};
    points     = varargin{2};
    dimension  = varargin{3};
    
else
    error('Too many input arguments.');
end
    
% shift zero-frequency to first bin 
shift            = zeros( 1, ndims(input_data) );
shift(dimension) = floor( size( input_data, dimension ) / 2 + 1 ); % -1

% buffer shift
input_data = circshift( input_data, shift );

% do FFT
out = fft( input_data, points, dimension );
%out = ifft( input_data, points, dimension, 'symmetric' );
%out = fft( real(input_data), points, dimension ) + ...
%    1j * fft( imag(input_data), points, dimension );

% normalization
out = out ./ (2*pi);
%out = (2*pi) .* out;

% buffer shift
out = circshift( out, shift );

end
