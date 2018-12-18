%% 初期設定（共通）
clear
dst_path = sprintf('H:/result/2018_12_06_IMCL_direct_estimation');
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/medium/2018_10_21_point_medium/point_medium_boundary_1.0mm_ IMCL0%.mat")
load("H:/data/kwave/result/2018_12_14_point_medium_various/boundary_1.0mm_IMCL0%/rfdata.mat")
load("H:/data/kwave/result/2018_12_14_point_medium_various/boundary_1.0mm_IMCL0%/kgrid.mat")
% 音速値
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
% IMCL割合
IMCL_rate = linspace(0,20,21);%[%]
num_IMCL = length(IMCL_rate);
v_muscle_with_IMCL = v_fat * IMCL_rate/100 + v_muscle*(1-IMCL_rate/100);%正解音速[m/s]
% 媒質境界位置
boundary_depth = linspace(19,0,20)*1e-3;
num_boundary_depth = length(boundary_depth);
ind_boundary_depth = zeros(num_boundary_depth,1);% kgrid上では境界位置はどのインデックスで表されるかをもとめる．
for i = 1:num_boundary_depth
    ind_boundary_depth(i) = find(single(kgrid.x_vec) == single(boundary_depth(i)));
end
num_boundary = length(ind_boundary_depth);
boundary_point = zeros(2,num_boundary);
boundary_point(1,:) = t_pos(1,51);
boundary_point(2,:) = kgrid.x_vec(ind_boundary_depth);
boundary_depth = (20 - linspace(19,0,20))*1e-3;% 実際の境界厚さ
% 探索位置
assumed_depth = 19e-3:-kgrid.dx:0;
num_assumed_depth = length(assumed_depth);
ind_assumed_depth = zeros(num_assumed_depth,1);% kgrid上では境界位置はどのインデックスで表されるかをもとめる．
for i = 1:num_assumed_depth
    ind_assumed_depth(i) = find(single(kgrid.x_vec) == single(assumed_depth(i)));
end
assumed_point = zeros(2,num_assumed_depth);
assumed_point(1,:) = t_pos(1,51);
assumed_point(2,:) = kgrid.x_vec(ind_assumed_depth);
% 素子配置
t_facing_distance = 0.04;%[m]
[~,num_receiver] = size(rfdata);
num_receiver = num_receiver/2;
element_pitch = abs(t_pos(1,1) - t_pos(1,2));
% 遅延プロファイル
distance_from_assumed_point = zeros(1,num_receiver);
% 推定値
assumed_SOS = v_muscle:-1:v_muscle_with_IMCL(end);
num_assumed_SOS = length(v_muscle_with_IMCL(end):v_muscle);
correlation = zeros(num_assumed_depth,num_assumed_SOS);
estimated_velocity = zeros(num_boundary,num_IMCL);
% 正解値
correct_velocity = zeros(num_boundary,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end
%% RFデータマスク処理部
loadpath = sprintf('H:/data/kwave/result/2018_12_13_layer_medium_various/boundary_5.0mm_IMCL10%%');
% loadpath = sprintf('H:/data/kwave/result/2018_12_14_point_medium_various/boundary_5.0mm_IMCL0%%');
load([loadpath,'/rfdata.mat'])
load([loadpath,'/kgrid.mat'])
load([loadpath,'/sourse_wave.mat'])
[num_sample,~] = size(rfdata);
% サンプル数を4倍にするためにスプライン補間
interp_sourcewave = interp1(source_wave,linspace(1,num_sample,num_sample*4)','spline');
interp_rfdata = zeros(num_sample*4,num_receiver);
rfdata_echo_only = zeros(num_sample*4,num_receiver);
delay_transmitted_wave = zeros(1,num_receiver);
% 透過波情報を消す
for jj = 1:num_receiver
    interp_rfdata(:,jj) = interp1(rfdata(:,jj),linspace(1,num_sample,num_sample*4),'spline');
    delay_transmitted_wave(1,jj) = round(((abs(assumed_point(1,1) - t_pos(1,jj)))/v_muscle)/(kgrid.dt/4));
    rfdata_echo_only(:,jj) = interp_rfdata(:,jj);
    rfdata_echo_only(1:delay_transmitted_wave(1,jj)+200,jj) = mean(rfdata_echo_only(1:delay_transmitted_wave(1,jj)+200,jj));
    rfdata_echo_only(:,jj) = abs(hilbert(rfdata_echo_only(:,jj)));
end
figure;
imagesc(rfdata_echo_only);
% 大津の二値化
BW_rfdata_echo_only = imbinarize(rfdata_echo_only);
figure;
imagesc(BW_rfdata_echo_only);
CC = bwconncomp(BW_rfdata_echo_only);
L = labelmatrix(CC);
S = regionprops(CC,'Area');
% 小さい検出領域を消す
mask_rfdata =  ismember(L, find([S.Area] >= max(struct2array(S))));
masked_rfdata = mask_rfdata .* interp_rfdata;
amp_masked_rfdata = abs(hilbert(masked_rfdata));
figure;
imagesc(amp_masked_rfdata);
% 送信波の包絡線
amp_sourcewave = abs(hilbert(interp_sourcewave));
auto_correlation_rfdata = diag(amp_masked_rfdata.' * amp_masked_rfdata);
auto_correlation_source_wave = amp_sourcewave.' * amp_sourcewave;
auto_correlation_source_wave = repmat(auto_correlation_source_wave,length(auto_correlation_rfdata),1);
delay_sourcewave = repmat(amp_sourcewave,1,num_receiver);
% 仮定遅延プロファイルと実測遅延プロファイルの相互相関を求める
figure;
for  kk = 1:num_assumed_depth
    for ll = 1:num_assumed_SOS
        for ii = 1:num_receiver
            distance_from_assumed_point(1,ii) = norm(assumed_point(:,kk)-t_pos(:,ii));%[m]
        end
        distance_round_trip = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
        delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
        for ii = 1:num_receiver
            source_wave2cat = amp_sourcewave;
            source_wave2cat(num_sample*4-delay_time_assumed(ii)+1:end,1)=NaN;
            source_wave2cat(isnan(source_wave2cat)) = [];
            delay_sourcewave(:,ii) = cat(1,zeros(delay_time_assumed(ii),1),source_wave2cat);
            correlation(kk,ll) = correlation(kk,ll) + (amp_masked_rfdata(:,ii).'*delay_sourcewave(:,ii)...
                /sqrt(auto_correlation_rfdata(ii,1)*auto_correlation_source_wave(ii,1)));
        end
        correlation(kk,ll) = correlation(kk,ll)/num_receiver;
    end
%     for ll = 1
%         for ii = 1:num_receiver
%             distance_from_assumed_point(1,ii) = norm(assumed_point(:,kk)-t_pos(:,ii));%[m]
%         end
%         distance_round_trip = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
%         delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
%         for ii = 1:num_receiver
%             source_wave2cat = amp_sourcewave;
%             source_wave2cat(num_sample*4-delay_time_assumed(ii)+1:end,1)=NaN;
%             source_wave2cat(isnan(source_wave2cat)) = [];
%             delay_sourcewave(:,ii) = cat(1,zeros(delay_time_assumed(ii),1),source_wave2cat);
%             correlation(kk,ll) = correlation(kk,ll) + (amp_masked_rfdata(:,ii).'*delay_sourcewave(:,ii)...
%                 /sqrt(auto_correlation_rfdata(ii,1)*auto_correlation_source_wave(ii,1)));
%         end
%         correlation(kk,ll) = correlation(kk,ll)/num_receiver;
%     end
%     imagesc(delay_sourcewave+amp_masked_rfdata);
%     caxis([0 0.1])
%     pause(0.2)
%     close gcf;
    disp(kk);
end
%% 保存部
% dst_path = sprintf('H:/result/2018_12_17_RF_data_mask_processing/2layer/2018_12_17_depth5.0mmIMCL0%%');
% if ~exist(dst_path, 'dir')
%     mkdir(dst_path);
% end
% figure;
% imagesc(IMCL_rate,20-assumed_depth*1e3,correlation);
% ylim([3 7])
% xlabel('IMCL content[%]')
% ylabel('depth[mm]')
% title('Correlation of delay curve')
% colorbar
% savefilename = sprintf('/correlation');
% savefig([dst_path,savefilename,'.fig'])
% exportfig([dst_path,savefilename],'png',[300,200])
