%%%%%%%%%%%%%%%%%%%%
% 対象：ワイヤ
% 波面遅延プロファイルにより音速推定
% 仮定遅延プロファイルと実測遅延プロファイルの相互相関
% を用いて音速を推定する
% 媒質   IMCL0%~20%(21種類)
% 深さ:
%%%%%%%%%%%%%%%%%%%%

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
boundary_point = zeros(2,num_boundary_depth);
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
distance_round_trip = zeros(1,num_receiver);
delay_time_assumed = zeros(1,num_receiver);
% 推定値
assumed_SOS = v_muscle:-1:v_muscle_with_IMCL(end);
num_assumed_SOS = length(v_muscle_with_IMCL(end):v_muscle);
correlation = zeros(num_assumed_depth,num_assumed_SOS);
estimated_velocity = zeros(num_boundary_depth,num_IMCL);
estimated_IMCL = zeros(num_boundary_depth,num_IMCL);
% 正解値
correct_velocity = zeros(num_boundary_depth,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end
%% 音速推定処理部
% 仮定遅延プロファイルと実測遅延プロファイルの相互相関を求める
for mm = 1:num_boundary_depth
    for nn = 1:num_IMCL
        loadpath = sprintf('H:/data/kwave/result/2018_12_14_point_medium_various/boundary_%0.1fmm_IMCL%d%%',...
            boundary_depth(mm)*1e3,IMCL_rate(nn));
        load([loadpath,'/rfdata.mat'])
        load([loadpath,'/kgrid.mat'])
        load([loadpath,'/sourse_wave.mat'])
        [num_sample,~] = size(rfdata);
        % rfデータ整形%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % サンプル数を4倍にするためにスプライン補間
        interp_sourcewave = interp1(source_wave,linspace(1,num_sample,num_sample*4)','spline');
        [~,offset_interp_sourcewave] = max(interp_sourcewave);% 送信波形の最大値を取る点．遅延曲線に散布図を重ね合わせることに使う．
        interp_rfdata = zeros(num_sample*4,num_receiver);
        rfdata_echo_only = zeros(num_sample*4,num_receiver);
        delay_transmitted_wave = zeros(1,num_receiver);
        for jj = 1:num_receiver
            interp_rfdata(:,jj) = interp1(rfdata(:,jj),linspace(1,num_sample,num_sample*4),'spline');
            delay_transmitted_wave(1,jj) = round(((abs(assumed_point(1,1) - t_pos(1,jj)))/v_muscle)/(kgrid.dt/4));
            rfdata_echo_only(:,jj) = interp_rfdata(:,jj);
            rfdata_echo_only(1:delay_transmitted_wave(1,jj)+200,jj) = mean(rfdata_echo_only(1:delay_transmitted_wave(1,jj)+200,jj));
        end
        auto_correlation_rfdata = diag(rfdata_echo_only.' * rfdata_echo_only);
        auto_correlation_source_wave = interp_sourcewave.' * interp_sourcewave;
        auto_correlation_source_wave = repmat(auto_correlation_source_wave,length(auto_correlation_rfdata),1);
        delay_sourcewave = repmat(interp_sourcewave,1,num_receiver);
        % 音速推定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for  kk = 1:num_assumed_depth
            for ll = 1:num_assumed_SOS
                for ii = 1:num_receiver
                    distance_from_assumed_point(1,ii) = norm(assumed_point(:,kk)-t_pos(:,ii));%[m]
                end
                distance_round_trip = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
                delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
                for ii = 1:num_receiver
                    source_wave2cat = interp_sourcewave;
                    source_wave2cat(num_sample*4-delay_time_assumed(ii)+1:end,1)=NaN;
                    source_wave2cat(isnan(source_wave2cat)) = [];
                    delay_sourcewave(:,ii) = cat(1,zeros(delay_time_assumed(ii),1),source_wave2cat);
                    correlation(kk,ll) = correlation(kk,ll) + (rfdata_echo_only(:,ii).'*delay_sourcewave(:,ii)...
                        /sqrt(auto_correlation_rfdata(ii,1)*auto_correlation_source_wave(ii,1)));
                end
                correlation(kk,ll) = correlation(kk,ll)/num_receiver;
            end
            dispname = sprintf('depth # is %d, IMCL # is %d',kk,nn);
            disp(dispname) %#ok<DSPS>
        end
        [~,ind_estimate_v] = max(max(correlation));
        [~,ind_estimate_d] = max(correlation(:,ind_estimate_v));
        estimated_velocity(mm,nn) = assumed_SOS(1,ind_estimate_v);
        estimated_IMCL(mm,nn) = 100*((v_muscle-estimated_velocity(mm,nn))/(v_muscle-v_fat));
        % delay sourcewave と rfdata を重ねて表示するための処理．%%%%%%
        for ii = 1:num_receiver
            distance_from_assumed_point(1,ii) = norm(assumed_point(:,ind_estimate_d)-t_pos(:,ii));%[m]
        end
        distance_round_trip = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
        delay_time_assumed = round(distance_round_trip / assumed_SOS(ind_estimate_v) / (kgrid.dt/4));%[sample]
        % 画像保存%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dst_path = sprintf('H:/result/2018_12_13_IMCL_direct_estimation_single_element/wire/2018_12_18/2018_12_18_depth%0.1fmmIMCL%d%%',...
            boundary_depth(mm)*1e3,IMCL_rate(nn));
        if ~exist(dst_path, 'dir')
            mkdir(dst_path);
        end
        figure;
        imagesc(IMCL_rate,20-assumed_depth*1e3,correlation);
        xlabel('IMCL content[%]')
        ylabel('depth[mm]')
        titlename = sprintf('IMCL: %d %%, depth: %0.1f mm',IMCL_rate(nn),boundary_depth(mm)*1e3);
        title({'Correlation of delay curve';titlename})
        colorbar
        savefilename = sprintf('/correlation');
        savefig([dst_path,savefilename,'.fig'])
        exportfig([dst_path,savefilename],'png',[300,200])
        close gcf
        figure;
        imagesc(rfdata_echo_only);
        hold on
        scatter(1:num_receiver,delay_time_assumed+offset_interp_sourcewave,'r.')
        caxis([0 0.2])
        ylim([min(delay_time_assumed) max(delay_time_assumed+2*offset_interp_sourcewave)])
        xlabel('element number')
        ylabel('time[sample]')
        axis square
        title(titlename)
        colorbar
        savefilename = sprintf('/solved_delay_curve');
        savefig([dst_path,savefilename,'.fig'])
        exportfig([dst_path,savefilename],'png',[350,300])
        close gcf
    end
end
dst_path2 = sprintf('H:/result/2018_12_13_IMCL_direct_estimation_single_element/wire/2018_12_18/figure');
if ~exist(dst_path2, 'dir')
    mkdir(dst_path2);
end
for mm = 1:num_boundary_depth
    figure;
    plot(correct_velocity(mm,:),estimated_velocity(mm,:));
    xlabel('correct velocity [m/s]')
    ylabel('estimated velocity [m/s]')
    titlename = sprintf('depth: %0.1f mm',boundary_depth(mm)*1e3);
    title(titlename)
    xlim([v_fat v_muscle])
    savefilename = sprintf('/velocity_depth_%0.1f mm',boundary_depth(mm)*1e3);
    savefig([dst_path2,savefilename,'.fig'])
    exportfig([dst_path2,savefilename],'png',[300,300])
    close gcf
    figure;
    plot(IMCL_rate(1,mm),estimated_IMCL(mm,:));
    xlabel('correct IMCL content [%]')
    ylabel('estimated IMCL content [%]')
    titlename = sprintf('depth: %0.1f mm',boundary_depth(mm)*1e3);
    title(titlename)
    xlim([0 20])
    ylim([0 20])
    savefilename = sprintf('/IMCL_depth_%0.1f mm',boundary_depth(mm)*1e3);
    savefig([dst_path2,savefilename,'.fig'])
    exportfig([dst_path2,savefilename],'png',[300,300])
    close gcf
end

dst_path3 = sprintf('H:/result/2018_12_13_IMCL_direct_estimation_single_element/wire/2018_12_18/');
if ~exist(dst_path3, 'dir')
    mkdir(dst_path3);
end
savefilename = sprintf('2018_12_18_all_result');
save([dst_path3,savefilename])