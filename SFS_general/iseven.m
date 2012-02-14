function bool = iseven(number)
%ISEVEN true for even integer
%
%   Usage: bool = iseven(number)
%
%   Input parameters:
%       number  - number (or vector/matrix with numbers) to be tested
%
%   Output parameters:
%       bool    - true if the number is even, false else
%
%   ISEVEN(number) checks if the given number is even.
%
% see also: isodd
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(number);


%% ===== Computation =====================================================

% Check if two is a factor of the given number
divided_by_two = mod( number,2 );

% Look for the even numbers
even_number_index = find( divided_by_two == 0 );

% Create answer
bool = false( size(number) );
bool(even_number_index) = true;
