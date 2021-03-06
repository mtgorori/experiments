%%%%%%%%%%%%%%
% 結果のfigureを作成する．2018/12/10
%%%%%%%%%%%%%%
load("H:\result\2018_12_06_IMCL_direct_estimation\2layer\2018_12_07_multi_layer_variable.mat")
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

% 全体図
figure;
plot(focal_depth*1e3,estimated_velocity);
hold on
plot(focal_depth*1e3,correct_velocity,'--');
hold off
xlabel('boundary depth[mm]');
ylabel('estimated velocity');
legend('estimated','correct');
savefilename = sprintf('/estimation_velocity');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% 4mm以降の拡大図
figure;
plot(focal_depth(3:end)*1e3,estimated_velocity(3:end));
hold on
plot(focal_depth(3:end)*1e3,correct_velocity(3:end),'--');
hold off
xlabel('boundary depth[mm]');
ylabel('estimated velocity');
legend('estimated','correct')
savefilename = sprintf('/estimation_velocity_detail');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% 深度ごとの推定誤差を棒グラフで示す
figure;
bar(focal_depth*1e3,estimated_velocity - correct_velocity);
xlabel('boundary depth[mm]');
ylabel('error of estimation[m/s] ');
savefilename = sprintf('/error_of_estimation');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% 4mm以降の拡大図
figure;
bar(focal_depth(3:end)*1e3,estimated_velocity(3:end) - correct_velocity(3:end));
xlabel('boundary depth[mm]');
ylabel('error of estimation[m/s] ');
savefilename = sprintf('/error_of_estimation_detail');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);