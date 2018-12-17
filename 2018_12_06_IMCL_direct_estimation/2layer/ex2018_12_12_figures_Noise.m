%%%%%%%%%%%%%%
% 結果のfigureを作成する．2018/12/12
%%%%%%%%%%%%%%
load("H:\result\2018_12_06_IMCL_direct_estimation\2layer\2018_12_07_multi_layer_variable.mat",...
    'estimated_velocity','aberration_ave_error','focal_depth','correct_velocity');
estimated_velocity_noNoise = estimated_velocity;
aberration_ave_error_noNoise = aberration_ave_error;
load("H:\result\2018_12_06_IMCL_direct_estimation\2layer\2018_12_12_Noise_signal_corr\20dB\2018_12_12_multi_layer_variable.mat",...
    'estimated_velocity','aberration_ave_error');
estimated_velocity_20dB = estimated_velocity;
aberration_ave_error_20dB = aberration_ave_error;
load("H:\result\2018_12_06_IMCL_direct_estimation\2layer\2018_12_12_Noise_signal_corr\10dB\2018_12_12_multi_layer_variable.mat",...
    'estimated_velocity','aberration_ave_error');
estimated_velocity_10dB = estimated_velocity;
aberration_ave_error_10dB = aberration_ave_error;
%%%%%%%%%%%%%%%%%%%%
% 対象：ワイヤターゲット，境界厚さ：2mm~19mm
% 設定音速：1580[m/s]＝正解音速
% 焦点水平位置固定：y=0
% 境界位置：既知
% IMCL割合を0 %に固定する．
% 波面遅延プロファイルにより音速推定
% 素子間受信波時間差は実信号の正規化相互相関を用いる
%%%%%%%%%%%%%%%%%%%%

%% 境界深さと音速推定精度の関係
dst_path = sprintf('H:/result/2018_12_06_IMCL_direct_estimation/2layer/2018_12_12_Noise_signal_corr');
% 全体図
figure;
plot(focal_depth*1e3,correct_velocity,'--');
hold on
plot(focal_depth*1e3,estimated_velocity_noNoise);
plot(focal_depth*1e3,estimated_velocity_20dB);
plot(focal_depth*1e3,estimated_velocity_10dB);
hold off
xlabel('boundary depth[mm]');
ylabel('estimated velocity');
legend('correct','no Noise','SNR=20 dB','SNR=10 dB');
savefilename = sprintf('/estimation_velocity');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% 4mm以降の拡大図
figure;
plot(focal_depth(3:end)*1e3,correct_velocity(3:end),'--');
hold on
plot(focal_depth(3:end)*1e3,estimated_velocity_noNoise(3:end));
plot(focal_depth(3:end)*1e3,estimated_velocity_20dB(3:end));
plot(focal_depth(3:end)*1e3,estimated_velocity_10dB(3:end));
hold off
xlabel('boundary depth[mm]');
ylabel('estimated velocity');
legend('correct','no Noise','SNR=20 dB','SNR=10 dB');
savefilename = sprintf('/estimation_velocity_detail');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
%RMS error
figure;
plot(focal_depth*1e3,sqrt(aberration_ave_error_noNoise)*1e9);
hold on
plot(focal_depth*1e3,sqrt(aberration_ave_error_20dB)*1e9);
plot(focal_depth*1e3,sqrt(aberration_ave_error_10dB)*1e9);
xlabel('boundary depth[mm]');
ylabel('RMS error[ns]');
ylim([80 500])
legend('no Noise','SNR=20 dB','SNR=10 dB','Location','northwest');
savefilename = sprintf('/ave_aberration_interp_rawsignal');
savefig([dst_path,savefilename,'.fig'])
exportfig([dst_path,savefilename],'png',[400,400])