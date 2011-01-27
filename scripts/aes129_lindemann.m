% Use the lindemann model to calculate binaural parameter of the stimuli
%
fs = 44100;
T = linspace(-1,1,45);
sprache = wavread('/home/hagen/d/data/signals/castanets.wav');

function [maxv,cent] = binaural_para(insig)
    T = linspace(-1,1,45);
    fs = 44100;
    [cc,t] = lindemann(insig(1:fs,:),fs,'stationary');
    [v,idx] = max(mean(cc,3));
    maxv = T(idx);
    tmp = mean(cc,3);
    cent = lindcentroid(tmp');
end

brir = wavread('brir_0m_xs0_ys1_R1_phi-60_fal2200.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L0_R1_60_max,L0_R1_60_c] = binaural_para(insig)
brir = wavread('brir_1m_xs0_ys1_R1_phi-60_fal2200.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L1_R1_60_max,L1_R1_60_c] = binaural_para(insig)
brir = wavread('brir_2m_xs0_ys1_R1_phi-60_fal2200.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L2_R1_60_max,L2_R1_60_c] = binaural_para(insig)
brir = wavread('brir_4m_xs0_ys1_R1_phi-60_fal2200.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L4_R1_60_max,L4_R1_60_c] = binaural_para(insig)

brir = wavread('brir_0m_xs0_ys1_R1_phi-30_fal4100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L0_R1_30_max,L0_R1_30_c] = binaural_para(insig)
brir = wavread('brir_1m_xs0_ys1_R1_phi-30_fal4100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L1_R1_30_max,L1_R1_30_c] = binaural_para(insig)
brir = wavread('brir_2m_xs0_ys1_R1_phi-30_fal4100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L2_R1_30_max,L2_R1_30_c] = binaural_para(insig)
brir = wavread('brir_2m_xs0_ys1_R4_phi-30_fal2100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L4_R1_30_max,L4_R1_30_c] = binaural_para(insig)

brir = wavread('brir_0m_xs0_ys1_R4_phi-60_fal1500.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L0_R4_60_max,L0_R4_60_c] = binaural_para(insig)
brir = wavread('brir_1m_xs0_ys1_R4_phi-60_fal1500.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L1_R4_60_max,L1_R4_60_c] = binaural_para(insig)
brir = wavread('brir_2m_xs0_ys1_R4_phi-60_fal1500.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L2_R4_60_max,L2_R4_60_c] = binaural_para(insig)
brir = wavread('brir_10m_xs0_ys1_R4_phi-60_fal1500.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L10_R4_60_max,L10_R4_60_c] = binaural_para(insig)

brir = wavread('brir_0m_xs0_ys1_R4_phi-30_fal2100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L0_R4_30_max,L0_R4_30_c] = binaural_para(insig)
brir = wavread('brir_1m_xs0_ys1_R4_phi-30_fal2100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L1_R4_30_max,L1_R4_30_c] = binaural_para(insig)
brir = wavread('brir_2m_xs0_ys1_R4_phi-30_fal2100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L2_R4_30_max,L2_R4_30_c] = binaural_para(insig)
brir = wavread('brir_10m_xs0_ys1_R4_phi-30_fal2100.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[L10_R4_30_max,L10_R4_30_c] = binaural_para(insig)

brir = wavread('brir_ref_xs0_ys1_R1_phi0.wav');
sigl = conv(sprache,brir(:,1));
sigr = conv(sprache,brir(:,2));
insig = [sigl sigr];
[ref_max,ref_c] = binaural_para(insig)


