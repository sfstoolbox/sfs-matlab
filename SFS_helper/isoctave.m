function bool=isoctave()
%ISOCTAVE  True if the operating environment is octave.
%   Usage: t=isoctave();
%
%   ISOCTAVE returns 1 if the operating environment is Octave, otherwise
%   0 (Matlab)

bool = ( exist('OCTAVE_VERSION','builtin') ~= 0 )
