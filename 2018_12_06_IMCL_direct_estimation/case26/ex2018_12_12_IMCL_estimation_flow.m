%%%%%%%%%%%%%%%%%%%%
% 対象：媒質case26
% 設定音速：1580[m/s]＝正解音速
% 焦点水平位置：y=0~4.8mm
% 境界位置：強度分布により推定
% IMCL割合を0 %に固定する．
% 波面遅延プロファイルにより音速推定
% 素子間受信波時間差は解析信号の正規化相互相関(FFTアプローチ)を用いる
%%%%%%%%%%%%%%%%%%%%

%% 初期設定（共通）
clear
dst_path = sprintf('H:/result/2018_12_06_IMCL_direct_estimation/case26/2018_12_12_flow');
if ~exist(dst_path, 'dir')
    mkdir(dst_path);
end
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/medium/2018_11_07_layer_medium/Layer_medium_boundary_2.0mm_ IMCL0%.mat")
load("H:/data/kwave/result/2018_11_07_layer_medium/Layer_medium_boundary_2.0mm_ IMCL0%/rfdata.mat")
load("H:/data/kwave/result/2018_11_07_layer_medium/Layer_medium_boundary_2.0mm_ IMCL0%/kgrid.mat")
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
v_reference = v_muscle;
rate_IMCL = linspace(1,20,20);
t_facing_distance = 0.04;%[m]
[~,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
num_rate_IMCL = length(rate_IMCL);
% num_medium  = num_rate_IMCL;
num_medium = 1;
reference_point = zeros(1,num_echo_receiver);
reference_point_lowerlimit = zeros(1,num_echo_receiver);%均質性評価のためのRFデータマスキングに使う
reference_point_upperlimit = zeros(1,num_echo_receiver);%均質性評価のためのRFデータマスキングに使う
distance_from_focal_point_all = zeros(1,num_echo_receiver);
delay_time_all = zeros(1,num_echo_receiver);
num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx/2 - 3;%'3'とあるのは，最近接距離が0.4 mmであることを考慮している．
focal_depth = zeros(1,num_depth);
focal_point = zeros(2,num_depth);
element_pitch = abs(t_pos(1,1) - t_pos(1,2));
lateral_range_max = 4.8*1e-3;%[m]
lateral_range_min = 0;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mmまで
ind_max_signal = zeros(num_lateral,num_medium);
max_signal = zeros(num_lateral,num_medium);
displacement = zeros(num_echo_receiver,num_lateral,num_medium);
error_wavefront = zeros(num_lateral,num_medium);
estimated_velocity = zeros(num_lateral,num_medium);
correct_velocity = v_muscle*(1-rate_IMCL/100)+v_fat*(rate_IMCL/100);
aberration_error = zeros(num_echo_receiver,num_lateral,num_medium);
aberration_ave_error = zeros(num_lateral,num_medium);
ind_min_error = zeros(1,num_medium);
min_error = zeros(1,num_medium);
focal_signal_total = zeros(num_depth,num_lateral,num_medium);
focal_signal = zeros(num_lateral,num_medium);
detected_boundary = zeros(num_lateral,num_medium);

%% 媒質の選択
for mm = 1:num_medium
    for ll = 1:num_lateral
        for ii = 9:num_depth
            focal_depth(1,ii) = single((ii+3)*kgrid.dx);
            focal_point(2,ii) = single(t_pos(2,1)-focal_depth(1,ii));
            focal_point(1,ii) = single(-(ll-1) * element_pitch);
        end
        pathname = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',mm);
        cd(pathname);
        load('rfdata.mat');
        load('kgrid.mat');
        load('sourse_wave.mat');
        [num_sample,~,~] = size(rfdata);
        [~,offset_source_wave] = max(abs(hilbert(source_wave)));
        %% 媒質境界を探索
        for ii = 9:num_depth
            %駆動素子の決定
            target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2)+focal_point(1,ii)));
            %受信用の参照点算出
            for jj = 1:num_echo_receiver
                distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                delay_time_all(1,jj) = floor(((distance_from_focal_point_all(1,jj) - focal_depth(1,ii))/v_reference)/kgrid.dt);%[sample]
                reference_point(1,jj) = floor(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference)/kgrid.dt+25-1);
            end
            %送信ビームフォーミング（共通）
            focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
            %受信ビームフォーカシング
            focal_signal_total(ii,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
        end
        [max_tmp,ind_tmp] = findpeaks(focal_signal_total(9:end,ll,mm),'Npeaks',1,'SortStr','descend');
        ind_max_signal(ll,mm) = ind_tmp+8;
        max_signal(ll,mm) = focal_signal_total(ind_max_signal(ll,mm),ll,mm);
        %% 検出した媒質境界にビームフォーカス
        %駆動素子の決定
        target_element = find(((-focal_depth(1,ind_max_signal(ll,mm))/2+focal_point(1,ind_max_signal(ll,mm)))...
            <=t_pos(1,1:100))&(t_pos(1,1:100)<=(focal_depth(1,ind_max_signal(ll,mm))/2+focal_point(1,ind_max_signal(ll,mm)))));
        %受信用の参照点算出
        for jj = 1:num_echo_receiver
            distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ind_max_signal(ll,mm)));
            delay_time_all(1,jj) = floor(((distance_from_focal_point_all(1,jj) - focal_depth(1,ind_max_signal(ll,mm)))/v_reference)/kgrid.dt);%[sample]
            reference_point(1,jj) = floor(delay_time_all(1,jj)+(2*focal_depth(1,ind_max_signal(ll,mm))/v_reference)/kgrid.dt+25-1);
        end
        %送信ビームフォーミング（共通）
        focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
        interp_rfdata = zeros(num_sample*4,num_echo_receiver);
        for jj = 1:num_echo_receiver
            interp_rfdata(:,jj) = interp(focused_rfdata(:,jj),4);
        end
        source_wave = interp(source_wave,4);
        %% 音速推定
        %RFデータマスキング（均質性評価のため)
        max_delay = (abs(diff(reference_point)) + 2)*4;
        min_reference = (min(reference_point) - 5)*4;
        max_reference = (min(reference_point) + 40 - 1)*4;
        mask_rfdata = zeros(size(interp_rfdata));
        mask_rfdata(min_reference:max_reference,:) = 1;
        focused_rfdata_mask = mask_rfdata .* interp_rfdata;
        % 波面形状推定
        delay_profile = zeros(1,length(target_element));
        for jj = 1:length(target_element)-1
            [acor,lag] = xcorr(interp_rfdata(:,target_element(jj+1)),focused_rfdata_mask(:,target_element(jj)),max_delay(target_element(jj)),'coeff');
            [~,I] = max(abs(acor));
            displacement(target_element(jj)+1,ll,mm) = lag(I);
            delay_profile(1,jj+1) = sum(displacement(target_element(1):target_element(jj)+1,ll,mm));
        end
        ind_central_element = find(target_element ==  round((min(target_element)+max(target_element))/2));
        delay_profile = delay_profile + abs(delay_profile(1,ind_central_element));
        [~,delay_offset] = findpeaks(abs(hilbert(focused_rfdata_mask(:,target_element(ind_central_element)))),'NPeaks',1,'SortStr','descend');
        [~,source_wave_offset] = findpeaks(abs(hilbert(source_wave)),'NPeaks',1,'SortStr','descend');
        delay_offset = delay_offset - source_wave_offset;
        delay_profile = delay_profile + delay_offset/2;
        poly_delay_profile_fitted = polyfit(t_pos(1,target_element),(delay_profile(1:length(target_element))*kgrid.dt/4).^2,2);
        estimated_velocity(ll,mm) = 1/sqrt(poly_delay_profile_fitted(1));
        f = polyval(poly_delay_profile_fitted.',t_pos(1,target_element).');
        T = table(t_pos(1,target_element).',(delay_profile*kgrid.dt/4).^2.',f,(delay_profile*kgrid.dt/4).^2.'-f,'VariableNames',{'X','Y','Fit','FitError'});
        aberration_error(target_element,ll,mm) = T.FitError;
        aberration_ave_error(ll,mm) = sum(abs(T.FitError))/length(target_element);
        %% 保存部
        dst_path2 = sprintf('H:/result/2018_12_06_IMCL_direct_estimation/case26/2018_12_12_flow/IMCL%d%%',mm);
        if ~exist(dst_path2, 'dir')
            mkdir(dst_path2);
        end
        figure;
        plot(focal_depth*1e3,focal_signal_total(:,ll,mm));
        hold on
        scatter(focal_depth(ind_max_signal(ll,mm))*1e3,max_signal(ll,mm),'red','filled')
        hold off
        xlabel('depth[mm]')
        ylabel('intensity[au]')
        savefilename1 = sprintf('/focal_signal_lateral%0.1fmm.fig',focal_point(1,9)*1e3);
        savefilename2 = sprintf('/focal_signal_lateral%0.1fmm',focal_point(1,9)*1e3);
        savefig([dst_path2,savefilename1])
        exportfig([dst_path2,savefilename2],'png',[400,400])
        close gcf
        delay_profile = zeros(1,length(target_element));
        for jj = 1:length(target_element)-1
            [acor,lag] = xcorr(interp_rfdata(:,target_element(jj+1)),focused_rfdata_mask(:,target_element(jj)),max_delay(target_element(jj)),'coeff');
            [~,I] = max(abs(acor));
            displacement(target_element(jj)+1,ll,mm) = lag(I);
            delay_profile(1,jj+1) = sum(displacement(target_element(1):target_element(jj)+1,ll,mm));
        end
        ind_central_element = find(target_element ==  round((min(target_element)+max(target_element))/2));
        delay_profile = delay_profile + abs(delay_profile(1,ind_central_element));
        [~,delay_offset] = findpeaks(abs(hilbert(focused_rfdata_mask(:,target_element(ind_central_element)))),'NPeaks',1,'SortStr','descend');
        [~,source_wave_offset] = findpeaks(abs(hilbert(source_wave)),'NPeaks',1,'SortStr','descend');
        delay_profile = delay_profile + delay_offset;
        figure;
        imagesc(t_pos(1,:)*1e3,kgrid.t_array*1e9,focused_rfdata_mask);
        hold on
        scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,delay_profile*kgrid.dt/4*1e9,'blue','filled');
        plot(t_pos(1,(min(target_element):max(target_element)))*1e3,sqrt(f),'red');
        hold off
        xlabel('lateral[mm]');
        ylabel('time[ns]');
        axis square;
        axis tight;
        xlim([t_pos(1,(min(target_element))-1)*1e3 t_pos(1,(max(target_element)+1))*1e3])
        ylim([min_reference*kgrid.dt*1e9/4 max_reference*kgrid.dt*1e9/4])
        colormap(bone);
        colorbar;
        savefilename = sprintf('/rfdata_detail_interp_rawsignal_lateral_%0.1fmm',focal_point(1,9)*1e3);
        savefig([dst_path2,savefilename,'.fig'])
        exportfig([dst_path2,savefilename],'png',[400,400])
        close gcf
        figure;
        imagesc(t_pos(1,:)*1e3,kgrid.t_array*1e9,abs(hilbert(interp_rfdata)));
        hold on
        scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,delay_profile*kgrid.dt/4*1e9,'blue','filled');
        % scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,(reference_point(1,(min(target_element):max(target_element)))+3)*kgrid.dt*1e9,'red');
        hold off
        xlabel('lateral[mm]');
        ylabel('time[ns]');
        axis square;
        axis tight;
        xlim([t_pos(1,(min(target_element))-1)*1e3 t_pos(1,(max(target_element)+1))*1e3])
        ylim([0 max_reference*kgrid.dt*1e9/2])
        colormap(bone);
        colorbar;
        caxis([0 max(max(abs(hilbert(focused_rfdata_mask))))])
        savefilename = sprintf('/rfdata_whole_interp_rawsignal_lateral_%0.1fmm',focal_point(1,9)*1e3);
        savefig([dst_path2,savefilename,'.fig'])
        exportfig([dst_path2,savefilename],'png',[400,400])
        close gcf
        figure;
        histogram(sqrt(abs(aberration_error(target_element,mm)))*1e9,'Normalization','probability','NumBins',round(sqrt(length(target_element))));
        ylim([0 0.8])
        xlabel('error of fitting[ns]');
        titlename = sprintf('lateral: %0.1f mm',focal_point(1,9)*1e3);
        title(titlename);
        savefilename = sprintf('/aberration_interp_rawsignal_lateral_%0.1fmm',focal_point(1,9)*1e3);
        savefig([dst_path2,savefilename,'.fig'])
        exportfig([dst_path2,savefilename],'png',[250,200])
        close gcf
    end
    loadpath = sprintf('H:/data/kwave/medium/2018_09_28_realisticScatter_variousIMCL/corrected/');
    loadfilename = sprintf('case26_IMCL%0.1f_pure.mat',mm);
    load([loadpath, loadfilename]);
    figure;
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    hold on
    scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
    plot(t_pos(2,1)*1e3 - focal_depth(ind_max_signal(:,mm))*1e3,...
        single(-(0:num_lateral-1) * element_pitch)*1e3,'r.-');
    axis equal
    axis tight
    xlabel('x-axis[mm]')
    ylabel('y-axis[mm]')
    caxis([1450 1580]);
    xlim([7.4 20.7])
    ylim([-9.1 4.2])
    c = colorbar;
    c.Label.String = '[m/s]';
    set(gca,'YDir','normal');
    savefilename = sprintf('/case26_IMCL%0.1f_pure_boundary',ll);
    savefigname = sprintf('/case26_IMCL%0.1f_pure_boundary.fig',ll);
    savefig([dst_path2,savefigname]);
    exportfig([dst_path2,savefilename],'png',[300 300]);
    close gcf
end
for mm = 1:num_medium
    [min_error(1,mm),ind_min_error(1,mm)] = min(aberration_ave_error(:,mm));
end
%% 保存部
figure;
scatter(1:num_medium,estimated_velocity(ind_min_error));
hold on
plot(correct_velocity,'--');
hold off
xlabel('IMCL content rarte[%]')
ylabel('estimated velocity[m/s]');
legend('estimated','correct')
savefilename = sprintf('/estimation_velocity');
savefig([dst_path,savefilename,'.fig'])
exportfig([dst_path,savefilename],'png',[600,500])
figure;
scatter(1:num_medium,sqrt(aberration_ave_error(ind_min_error))*1e9);
xlabel('IMCL content rarte[%]');
ylabel('RMS error[ns]');
savefilename = sprintf('/ave_aberration_interp_rawsignal');
savefig([dst_path,savefilename,'.fig'])
exportfig([dst_path,savefilename],'png',[300,250])
save([dst_path,'\2018_12_07_multi_layer_variable']);