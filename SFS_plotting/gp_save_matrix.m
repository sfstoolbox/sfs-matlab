function gp_save_matrix(file,x,y,M)
% GP_SAVE_MATRIX save x,y,M in the matrix binary format of Gnuplot
%   Usage: gp_save_matrix(file,x,y,M)
%
%   Input parameters:
%       file    - filename of the data file
%       x       - x axis values
%       y       - y axis values
%       M       - matrix data size(M)=x,y
%
%   GP_SAVE_MATRIX(file,x,y,M) saves the values of x,y and M in a binary matrix
%   format useable by Gnuplot, see
%   http://www.gnuplot.info/docs_4.4/gnuplot.html#x1-33700077.1.1 for details.
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
error(nargchk(4,4,nargin));

if ~ischar(file)
    error('%s: file has to be a string!',upper(mfilename));
end
if ~isnumeric(x) || ~isvector(x)
    error('%s: x has to be a vector!',upper(mfilename));
end
if ~isnumeric(y) || ~isvector(y)
    error('%s: y has to be a vector!',upper(mfilename));
end
if ~isnumeric(M) || ~ismatrix(M)
    error('%s: M has to be a matrix!',upper(mfilename));
end


%% ===== Computation =====================================================

% Check if the data has the right format
[lx,ly] = size(M);
if lx~=length(x) || ly~=length(y)
    error('%s: size(M) has to be x,y!',upper(mfilename));
end

% Create matrix to store in the file
MS = zeros(length(x)+1,length(y)+1);
MS(1,1) = length(y);
MS(1,2:end) = y;
MS(2:end,1) = x;
MS(2:end,2:end) = M;

% Write data into the file
fid = fopen(file,'w');
fwrite(fid,MS,'float');
fclose(fid);
