function [ out ] = fftx( varargin )
%FFTX Spatial fourier transform
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

% shift zero-frequency to first bin as ifft expects it
shift            = zeros( 1, ndims(input_data) );
shift(dimension) = floor( size( input_data, dimension ) / 2 - 1 );

% buffer shift
input_data = circshift( input_data, shift );

% do FFT
out = ifft( input_data, points, dimension );
%out = ifft( real(input_data), points, dimension, 'symmetric' ) + ...
%    1j * ifft( imag(input_data), points, dimension, 'symmetric' );

% normalization
out = 2*pi .* out;
%out = out ./ size(out, dimension);

% buffer shift
out = circshift( out, shift );

    
end
