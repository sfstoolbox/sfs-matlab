function status = test_all(modus)
%TEST_ALL runs all tests
%
%   Usage: test_all(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false

%   TEST_ALL(modus) runs all test function that should be checked before
%   preparing a new release of the toolbox.

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


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Main ============================================================
if modus
    disp('Running graphical modus, please inspect the figures');
    disp('');
    disp('Running "test_colormaps(1)"');
    test_colormaps(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_delayline(1)"');
    test_delayline(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_driving_functions_imp_with_delay(1)"');
    test_driving_functions_imp_with_delay(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_driving_functions(1)"');
    test_driving_functions(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_hrtf_extrapolation(1,''QU_KEMAR'')"');
    test_hrtf_extrapolation(1,'QU_KEMAR');
    disp('Hit Enter to continue');
    pause
    disp('Running "test_interpolation_methods(1)"');
    test_interpolation_methods(1);
    pause
    disp('Running "test_interpolation_point_selection(1)"');
    test_interpolation_point_selection(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_localwfs_vss(1)"');
    test_localwfs_vss(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_non_regular_grid(1)"');
    test_non_regular_grid(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_plot(1)"');
    test_plotting(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_secondary_source_positions(1)"');
    test_secondary_source_positions(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_wfs_25d(1)"');
    test_wfs_25d(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_imp_25d(1)"');
    test_imp_25d(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_fft_ifft(1)"');
    test_spectrum_signal_conversion(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_modal_weighting(1)"');
    test_modal_weighting(1);
    disp('Hit Enter to continue');
    pause
    disp('Running "test_sphbesselh_zeros(1)"');
    test_sphbesselh_zeros(1);
else
    if ~all([test_binaural_synthesis(0); ...
             test_delayline(0); ...
             test_driving_functions_imp_with_delay(0); ...
             test_driving_functions(0); ...
             test_hrtf_extrapolation(0,'QU_KEMAR'); ...
             test_interpolation_point_selection(0); ...
             test_localwfs_vss(0); ...
             test_non_regular_grid(0); ...
             test_secondary_source_diameter(0); ...
             test_secondary_source_positions(0); ...
             test_secondary_source_selection(0); ...
             test_tapering_window(0); ...
             test_wfs_25d(0); ...
             test_imp_25d(0); ...
             test_spectrum_signal_conversion(0);
             test_modal_weighting(0);
             test_sphbesselh_zeros(0);
             %test_wfs_iir_prefilter(0); ... % needs DSP Tooblox
            ])
        return;
    end
end


status = true;
