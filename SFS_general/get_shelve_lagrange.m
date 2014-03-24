function [H] = get_shelve_lagrange(f,H,FlagSub,fSub,FlagAliasing,fAliasing,Bandwidth_in_Oct)
%
%   this function will be called from hpre = wfs_iir_prefilter(conf)
%
%   note 1: there is a analytical expression for arbitrary Bandwidth_in_Oct
%   which is not implemented yet
%   note 2: the Offset_dB specified for Bandwidth_in_Oct work only for the
%   3dB/oct. slope, other driving functions may have other slopes and
%   different interpolation offsets may be required
%
%   References:
%       F. Schultz, V. Erbes, S. Spors, S. Weinzierl (2013) - "Derivation
%       of IIR prefilters for soundfield synthesis using linear secondary
%       source distributions", In: Proc. of the International Conference
%       on Acoustics AIA-DAGA, p.2372-2375

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
        %aliasing frequency limitation:
        tmp=(find(f<=fAliasing)); falias_i = tmp(end);
        H_Alias = H*0+H(falias_i);
        H_Alias(1) = H_Alias(2);
    end
    
    if FlagSub
        %Lagrange interpolation for SUB:
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

        %Lagrange interpolation tmp variables:
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
        %Lagrange interpolation for ALIASING:
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

        %Lagrange interpolation tmp variables:
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
