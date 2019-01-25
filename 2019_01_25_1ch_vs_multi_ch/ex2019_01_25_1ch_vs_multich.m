load("H:\result\2019_01_25_1ch_vs_multi_ch\v_and_v_reference.mat")

v_1ch = v(50,:,:);
v_1ch = reshape(v_1ch,model_num,freq_num);
v_1ch = v_1ch';
estimated_IMAT_1ch = ((v_1ch-1580)/(1450-1580))*100;

true_IMAT_1ch =  v_reference(50,:,:);
true_IMAT_1ch = reshape(true_IMAT_1ch,model_num,freq_num);
true_IMAT_1ch = true_IMAT_1ch';
true_IMAT_1ch = ((true_IMAT_1ch-1580)/(1450-1580))*100;

cd 'H:\data\kwave\result\2018_10_15_variousFrequency_Correct'
csvwrite('true_IMAT_1ch.csv',true_IMAT_1ch);
csvwrite('estimated_IMAT_1ch.csv',estimated_IMAT_1ch);