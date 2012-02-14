function bool = isodd(number)
%ISODD true for odd integer
%
%   Usage: bool = isodd(number)
%
%   Input parameters:
%       number  - number (or vector/matrix with numbers) to be tested
%
%   Output parameters:
%       bool    - true if the number is odd, false else
%
%   ISODD(number) checks if the given number is odd.
%
% see also: iseven
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(number)


%% ===== Computation =====================================================

% Check if two is a factor of the given number
divided_by_two = mod( number,2 );

% Look for the even numbers
odd_number_index = find( divided_by_two == 1 );

% Create answer
bool = false( size(number) );
bool(odd_number_index) = true;
