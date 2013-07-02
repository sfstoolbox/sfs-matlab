clear all; close all; clc;
%see wfs_iir_prefilter.m for details
SFS_start;
SFS_version
conf.fs = 44100;
conf.hpreflow = 200;
conf.hprefhigh = 1500;
conf.hpreBandwidth_in_Oct = 2;  
conf.hpreIIRorder = 4; 
hpre1 = wfs_iir_prefilter()      %call with no conf to test default values
hpre2 = wfs_iir_prefilter(conf)  %call with parameters in conf
SFS_stop;