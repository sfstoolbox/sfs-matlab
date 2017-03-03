function H = get_shelve_lagrange(f,H,FlagSub,fSub,FlagAliasing,fAliasing,Bandwidth_in_Oct)
%GET_SHELVE_LAGRANGE does an Lagrange interpolation towards shelving filter
%
%   Usage: H = get_shelve_lagrange(f,H,FlagSub,fSub,FlagAliasing, ...
%                                  fAliasing,Bandwidth_in_Oct)
%
%   Input parameters:
%       f                -  frequency vector in Hz, typical 0 to half
%                           sampling frequency, equidistant sampling is
%                           assumed (i.e. DFT frequencies) but not required
%       H                -  complex spectrum (mag/phase) at specified frequencies
%                           with meaningful slope, typically used for +3dB/oct.
%       FlagSub          -  use low-shelf part, 0 or 1
%       fSub             -  cut-frequency in Hz for low-shelf,
%                           check e.g. fLow=100Hz in [Fig2, Sch13]
%       FlagAliasing     -  use high-sehlf part, 0 or 1
%       fAliasing        -  cut-frequency in Hz for high-shelf,
%                           check e.g. fAliasing=2kHz in [Fig2, Sch13]
%       Bandwidth_in_Oct -  interpolation bandwidth, i.e. shelf knee range
%                           in octaves, allowed: 0.5, 1, 2, 3, 4
%
%   Output parameters:
%       H                -  return complex spectrum at specified frequencies
%                           after low-/ and high-shelf interpolation              -
%
%   GET_SHELVE_LAGRANGE(f,H,FlagSub,fSub,FlagAliasing,fAliasing,Bandwidth_in_Oct)
%   does an Lagrange interplation to get a shelving filter. This function is used
%   in wfs_iir_prefilter(conf).
%
%   Note 1: there is a analytical expression for arbitrary Bandwidth_in_Oct
%   which is not implemented yet
%   Note 2: the Offset_dB specified for Bandwidth_in_Oct work only for the
%   3dB/oct. slope, other driving functions may have other slopes and
%   different interpolation offsets may be required
%
%   References:
%       [Sch13]   F. Schultz, V. Erbes, S. Spors, S. Weinzierl (2013) - 
%       "Derivation of IIR prefilters for soundfield synthesis using linear
%       secondary source distributions", In: Proc. of the International
%       Conference on Acoustics AIA-DAGA, p.2372-2375
%
%   See also: wfs_iir_prefilter

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
% Revision: 07/02/2013 frank.schultz@uni-rostock.de initial development      *
%*****************************************************************************
    Hphase = unwrap(angle(H));  %save original phase
    H = 20*log10(abs(H));       %interpolation in dB

    if Bandwidth_in_Oct==4
        Offset_dB = 1.5053;
    elseif Bandwidth_in_Oct==3
        Offset_dB = 1.1286;
    elseif Bandwidth_in_Oct==2
        Offset_dB = 0.7528;
    elseif Bandwidth_in_Oct==1
        Offset_dB = 0.3766;
    elseif Bandwidth_in_Oct==0.5
        Offset_dB = 0.1857;
    else
        disp('error, Bandwidth_in_Oct must be 1,2,3 or 4')
    end

    if FlagSub
        %Sub-Bass limitation:
        tmp=(find(f<=fSub)); fSub_i = tmp(end);
        H_Sub = H*0+H(fSub_i);
        H_Sub(1) = H_Sub(2);
    end

    if FlagAliasing
        % Aliasing frequency limitation:
        tmp=(find(f<=fAliasing)); falias_i = tmp(end);
        H_Alias = H*0+H(falias_i);
        H_Alias(1) = H_Alias(2);
    end

    if FlagSub
        % Lagrange interpolation for SUB:
        fl = f(fSub_i)*2^(-Bandwidth_in_Oct/2);
        fh = f(fSub_i)*2^(+Bandwidth_in_Oct/2);
        tmp=(find(f<=fl)); fl_i = tmp(end);
        tmp=(find(f<=fh)); fh_i = tmp(end);

        P1x = log10(f(fl_i));
        P1y = H_Sub(fl_i);
        P2x = log10(f(fSub_i));
        P2y = H_Sub(fSub_i)+Offset_dB;
        P3x = log10(f(fh_i));
        P3y = H(fh_i);

        % Lagrange interpolation tmp variables:
        ipol1=P1y/((P1x-P2x)*(P1x-P3x));
        ipol2=P2y/((P2x-P1x)*(P2x-P3x));
        ipol3=P3y/((P3x-P1x)*(P3x-P2x));

        for fi=1:fl_i-1
        H(fi) = H_Sub(fi);
        end

        for fi=fl_i:fh_i
            H(fi)=          (ipol1*(log10(f(fi))-P2x)*(log10(f(fi))-P3x))+...
                            (ipol2*(log10(f(fi))-P1x)*(log10(f(fi))-P3x))+...
                            (ipol3*(log10(f(fi))-P1x)*(log10(f(fi))-P2x));
        end
    end

    if FlagAliasing
        % Lagrange interpolation for ALIASING:
        fl = f(falias_i)*2^(-Bandwidth_in_Oct/2);
        fh = f(falias_i)*2^(+Bandwidth_in_Oct/2);
        tmp=(find(f<=fl)); fl_i = tmp(end);
        tmp=(find(f<=fh)); fh_i = tmp(end);

        P1x = log10(f(fl_i));
        P1y = H(fl_i);
        P2x = log10(f(falias_i));
        P2y = H_Alias(falias_i)-Offset_dB;
        P3x = log10(f(fh_i));
        P3y = H_Alias(fh_i);

        % Lagrange interpolation tmp variables:
        ipol1=P1y/((P1x-P2x)*(P1x-P3x));
        ipol2=P2y/((P2x-P1x)*(P2x-P3x));
        ipol3=P3y/((P3x-P1x)*(P3x-P2x));

        for fi=fl_i:fh_i
            H(fi)=      (ipol1*(log10(f(fi))-P2x)*(log10(f(fi))-P3x))+...
                        (ipol2*(log10(f(fi))-P1x)*(log10(f(fi))-P3x))+...
                        (ipol3*(log10(f(fi))-P1x)*(log10(f(fi))-P2x));
        end

        for fi=fh_i+1:length(H)
            H(fi) = H_Alias(fi);
        end
    end

    %apply original phase and delog
    H = 10.^(H/20).*exp(1i*Hphase);
    H(1) = abs(H(2));
    %we ignore the special treatment for DC and fs/2 as would be usual for
    %DFT data, since we need a 'classical' frequency response
end
