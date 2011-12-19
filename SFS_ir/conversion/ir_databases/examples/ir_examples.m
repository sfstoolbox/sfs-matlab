function ir_examples()

i = -90;
while i <= 90
	irs = read_irs('QU_KEMAR_Auditorium3_1_src1.mat'); 
	ir = get_ir(irs,correct_azimuth(rad(i)));
	sig = auralize_ir(ir,'speech');
	j = i
	if i < 0
		j = i +360;
	end
	wavfile = sprintf('AUDITORIUM_1_1_%d.wav',j);
	wavwrite (sig, 44100, wavfile);
	i= i + 45 ;
end

i = -142;
while i <= 38
	irs = read_irs('QU_KEMAR_Auditorium3_1_src2.mat'); 
	ir = get_ir(irs,correct_azimuth(rad(i)));
	sig = auralize_ir(ir,'speech');
	j = i
	if i < 0
		j = i +360;
	end
	wavfile = sprintf('AUDITORIUM_1_2_%d.wav',j);
	wavwrite (sig, 44100, wavfile);
	i= i + 45 ;
end

i = 139;
while i <= 319
	irs = read_irs('QU_KEMAR_Auditorium3_1_src3.mat'); 
	ir = get_ir(irs,correct_azimuth(rad(i)));
	sig = auralize_ir(ir,'speech');
	j = i
	if i < 0
		j = i +360;
	end
	wavfile = sprintf('AUDITORIUM_1_3_%d.wav',j);
	wavwrite (sig, 44100, wavfile);
	i= i + 45 ;
end






i = -90;
while i <= 90
	irs = read_irs('QU_KEMAR_Auditorium3_2_src1.mat'); 
	ir = get_ir(irs,correct_azimuth(rad(i)));
	sig = auralize_ir(ir,'speech');
	j = i
	if i < 0
		j = i +360;
	end
	wavfile = sprintf('AUDITORIUM_2_1_%d.wav',j);
	wavwrite (sig, 44100, wavfile);
	i= i + 45 ;
end

i = -60;
while i <= 120
	irs = read_irs('QU_KEMAR_Auditorium3_2_src2.mat'); 
	ir = get_ir(irs,correct_azimuth(rad(i)));
	sig = auralize_ir(ir,'speech');
	j = i
	if i < 0
		j = i +360;
	end
	wavfile = sprintf('AUDITORIUM_2_2_%d.wav',j);
	wavwrite (sig, 44100, wavfile);
	i= i + 45 ;
end

i = -120;
while i <= 60
	irs = read_irs('QU_KEMAR_Auditorium3_2_src3.mat'); 
	ir = get_ir(irs,correct_azimuth(rad(i)));
	sig = auralize_ir(ir,'speech');
	j = i
	if i < 0
		j = i +360;
	end
	wavfile = sprintf('AUDITORIUM_2_3_%d.wav',j);
	wavwrite (sig, 44100, wavfile);
	i= i + 45 ;
end
