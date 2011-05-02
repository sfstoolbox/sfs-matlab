

Implementation notes:

=====================================================================================================

In the folder NFC_HOA_Filter_Implementation are all desired MATLAB function to achieve the following
filter designs for near-field compensated higher order ambisonics:

(1) FIR filter approximated by fir2.m method from MATLAB

(2) FIR filter approximatd by firintpolc.m method from Hans W. Schüssler

(3) IIR filter approximated by bilinear transformation

=====================================================================================================

The four steps to get the impulse responses of the FIR Filter (1) or (2):

- choose the calculation method of the driving function (with recursive or non recursive calculation 
  of the hankel function)

- set the right parameter at the section "filter design" (for example "E = pw_b_r;")

- choose the desired delay line at the section "desired delay line for filter"

- then choose the desired design method (1) or (2)


For the IIR filter design you only have to choose the IIR filter design (3) for plane wave or point
source.

=====================================================================================================


If you want to get a various filter design -for example for different incidence angle of a plane wave-
you have to change the parameters in HOA_main.m and do the steps mentioned above.
For example config.nLS = 28, config.fs = 96000 , ... .

The impulse responses will be stored in the folder IRs which is a subfolder of 
NFC_HOA_Filter_Implementation.
If you want to change properties like subfolder or name of the achieved impulse responses you have
to set your own properties in HOA_main.m (for example wavwrite('your impulse responses', 'sampling frequency,'your_folder'\'your_name.wav');).

=====================================================================================================

If you want to plot the reproduced wavefield in time or frequency domain open 

"plot_wave_field_from_imp_resp.m" to plot the reproduced wavefield in frequency domain

			or

"reproduced_field_time_domain_IR.m" to plot reproduced wavefield in time domain

Important note: 
You have to change the name of the *.wav-file (and/or the name of the subfolder) in
the "wavread()" method in the previous mentioned functions in order to get the wavefield-plot that 
you like.

=====================================================================================================

