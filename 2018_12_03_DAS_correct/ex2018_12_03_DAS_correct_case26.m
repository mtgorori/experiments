%%%%%%%%%%%%%%%%%%%%%%%%%%%%％％％％％％％%%%%%
% 対象：case26
% 焦点水平位置固定：y=0
% 境界面を検出→参照点の妥当性評価
% ch 50直下のピクセルにフォーカスをかける．
% 参照点：整相加算における参照点
% 比較する点：各chでのrf信号最大値(絶対値はとらない)
%誤検出が見られるが，反射強度を考慮して誤検出と判別できるようにする．
% フォーカス深度を１グリッドずつ変化させて，開口合成受信データを生成する．
% F値を固定する．最近接距離をいくつに設定する？
% x-axis: 20 mmのときに全素子を使うという前提を設ける．
% 対象にする媒質データ：case26を用いる．
% IMCL割合を0 %に固定する．
% kgrid.x_vec が0となるのは251番目の要素．
% グリッド幅が0.1 mm
% 素子間ピッチが0.4 mm
% よって，最近接深度(最近接距離は0.4 mm)
% 深度ごとの細かい素子割当は，floor()を用いる．
% 整相加算の際に参照する音速を1580 m/sとする．
% 整相加算の前に参照点と受信chごとの振幅(ヒルベルト変換後絶対値)最大値との距離の
% 合計を評価関数として媒質の均質性を評価することも同時に行う．
% update:rf-dataの配列サイズが大きくなってきたので，変数区分を細分化して各変数の呼び出し速度を上げる．[2018/11/05]
% update:for文ごとにフォルダを作成するようにする．フォルダごとにrfデータを保存する．
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 初期設定に必要なデータの呼び出し
clear;
load("H:/data/kwave/config/t_pos_2board.mat");
pathname = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',1.0);
cd(pathname);
load('rfdata.mat');
load('kgrid.mat');

%% 初期設定
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
rate_IMCL = 1;
v_reference = zeros(1,length(rate_IMCL));
t_facing_distance = 0.04;%[m]
[~,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
num_rate_IMCL = length(rate_IMCL);
num_medium  = num_rate_IMCL;
reference_point = zeros(1,num_echo_receiver);
reference_point_lowerlimit = zeros(1,num_echo_receiver);%均質性評価のためのRFデータマスキングに使う
reference_point_upperlimit = zeros(1,num_echo_receiver);%均質性評価のためのRFデータマスキングに使う
point_max_in_mask = zeros(1,num_echo_receiver);%マスク処理後のRFデータで振幅最大のサンプル点情報
distance_from_focal_point_all = zeros(1,num_echo_receiver);

num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx/2 - 3;%'3'とあるのは，最近接距離が0.4 mmであることを考慮している．
focal_depth = zeros(1,num_depth);
focal_point = zeros(2,num_depth);
element_pitch = abs(t_pos(1,1) - t_pos(1,2));
lateral_range_max = 10.0*1e-3;
lateral_range_min = 0;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mmまで

ind_max_signal = zeros(num_rate_IMCL,num_medium,num_lateral);
max_signal = zeros(num_rate_IMCL,num_medium,num_lateral);
ind_search_depth = zeros(num_lateral,num_medium);%探索範囲を一回目の推定で決める．その際の探索範囲の中央値
displacement = zeros(num_echo_receiver,num_rate_IMCL,num_medium,num_lateral);
displacement_ref = zeros(num_echo_receiver,num_rate_IMCL,num_medium,num_lateral);
error_wavefront = zeros(num_rate_IMCL,num_medium,num_lateral);
focal_signal_total = zeros(num_depth,num_rate_IMCL,num_medium,num_lateral);
focal_signal = zeros(num_rate_IMCL,num_medium,num_lateral);
detected_boundary = zeros(num_rate_IMCL,num_medium,num_lateral);

%% 処理部
for mm = 1:num_lateral
    for ii = 9:num_depth
        focal_depth(1,ii) = single((ii+3)*kgrid.dx);
        focal_point(2,ii) = single(t_pos(2,1)-focal_depth(1,ii));
        focal_point(1,ii) = single((mm-1) * element_pitch);
    end
    for ll = 1:num_medium
        %% 参照配列の呼び出し
        pathname = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',ll);
        cd(pathname);
        load('rfdata.mat');
        load('kgrid.mat');
        [num_sample,~,~] = size(rfdata);
        %% 境界位置の候補を選定する．(候補数：num_IMCL x num_medium)
        for kk = 1:num_rate_IMCL
            %% 境界位置のアタリをつける(full depth)
            if kk == 1
                %% 反射強度プロファイルを求める．
                v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
                for ii = 9:num_depth
                    target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                    %受信用の参照点算出
                    for jj = 1:num_echo_receiver
                        distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                        delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                        reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    end
                    %送信ビームフォーミング（共通）
                    focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                    %受信ビームフォーカシング
                    focal_signal_total(ii,kk,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
                end
                %% 反射強度が最大の焦点位置を探索．
                [max_tmp,ind_tmp] = findpeaks(abs(focal_signal_total(9:end,kk,ll,mm)),'Npeaks',1,'SortStr','descend');
                [min_ind_tmp,ind_ind_tmp] = min(ind_tmp);
                ind_max_signal(kk,ll,mm) = min_ind_tmp+8;
                max_signal(kk,ll,mm) = focal_signal_total(ind_max_signal(kk,ll,mm),kk,ll,mm);
                ind_search_depth(mm,ll) = ind_max_signal(kk,ll,mm);
                ii = ind_search_depth(mm,ll);
                %% 参照点と波面との形状相関を求める．（相互相関法）
                target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                %受信用の参照点算出
                for jj = 1:num_echo_receiver
                    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                    delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                    reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    %25はfocal_amplitudeを最大にするオフセット．
                end
                max_delay = abs(diff(reference_point)) + 2;
                %送信ビームフォーミング（共通）
                focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                %RFデータマスキング（均質性評価のため)
                min_reference = min(reference_point) - 5;
                max_reference = min(reference_point) + 40 - 1;
                mask_rfdata = zeros(size(focused_rfdata));
                mask_rfdata(min_reference:max_reference,:) = 1;
                focused_rfdata_masked = mask_rfdata .* focused_rfdata;
                for jj = 1:length(target_element)-1
                    %波面形状評価
                    displacement(target_element(jj),kk,ll,mm)  = - finddelay(focused_rfdata(:,target_element(jj+1)),...
                        focused_rfdata_masked(:,target_element(jj)),max_delay(target_element(jj)));
                end
                displacement_ref(target_element(1:end-1),kk,ll,mm) = diff(reference_point(target_element));
                error_wavefront(kk,ll,mm) = sum(abs(displacement(target_element(1:end-1),kk,ll,mm) - ...
                    displacement_ref(target_element(1:end-1),kk,ll,mm)))/(length(target_element)-1);
                %% 出力データの保存
                dst_path = sprintf('H:/result/2018_12_03_DAS_correct/case26/correct/lateral%0.1fmm/true%d_assumption%d',...
                    focal_point(1,9)*1e3,rate_IMCL(1,ll),rate_IMCL(1,kk));
                if ~exist(dst_path, 'dir')
                    mkdir(dst_path);
                end
                displacement_tmp = displacement(:,kk,ll,mm);
                save([dst_path,'\pre_rfdata.mat'],'focused_rfdata','focused_rfdata_masked','displacement_tmp','reference_point');
                figure;
                imagesc(focused_rfdata_masked);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),'red');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                xlim([min(target_element) max(target_element)])
                ylim([min_reference max_reference])
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_masked.fig'])
                exportfig([dst_path,'\pre_rfdata_masked'],'png',[400,400])
                close gcf
                figure;
                imagesc(focused_rfdata);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),2,'red');
                scatter(min(target_element):max(target_element),point_max_in_mask(min(target_element):max(target_element)),2,'blue','filled');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_whole.fig'])
                exportfig([dst_path,'\pre_rfdata_whole'],'png',[400,400])
                close gcf
                figure;
                plot(focal_depth*1e3,focal_signal_total(:,kk,ll,mm));
                hold on
                scatter(focal_depth(ind_max_signal(kk,ll,mm))*1e3,max_signal(kk,ll,mm),'red','filled')
                hold off
                xlabel('焦点深さ[mm]')
                ylabel('信号強度[au]')
                savefig([dst_path,'\pre_focal_signal.fig'])
                exportfig([dst_path,'\pre_focal_signal'],'png',[400,400])
                close gcf
                fprintf('prepareing... lateral number is %d, medium number is %d, and estimation number is %d\n',mm,ll,kk);
            else
            %% 候補を絞って境界位置を求める．(limited depth)
                %% 反射強度プロファイルを求める．
                v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
                for ii = ind_search_depth(mm,ll)-5:ind_search_depth(mm,ll)+5
                    target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                    %受信用の参照点算出
                    for jj = 1:num_echo_receiver
                        distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                        delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                        reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    end
                    %送信ビームフォーミング（共通）
                    focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                    %受信ビームフォーカシング
                    focal_signal_total(ii,kk,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
                end
                %% 反射強度が最大の焦点位置を探索．
                [max_tmp,ind_tmp] = findpeaks(abs(focal_signal_total(9:end,kk,ll,mm)),'Npeaks',1,'SortStr','descend');
                [min_ind_tmp,ind_ind_tmp] = min(ind_tmp);
                ind_max_signal(kk,ll,mm) = min_ind_tmp+8;
                max_signal(kk,ll,mm) = focal_signal_total(ind_max_signal(kk,ll,mm),kk,ll,mm);
                ii = ind_max_signal(kk,ll,mm);
                %% 参照点と波面との形状相関を求める．（相互相関法）
                target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                %受信用の参照点算出
                for jj = 1:num_echo_receiver
                    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                    delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                    reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    %25はfocal_amplitudeを最大にするオフセット．
                end
                max_delay = abs(diff(reference_point)) + 2;
                %送信ビームフォーミング（共通）
                focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                %RFデータマスキング（均質性評価のため)
                min_reference = min(reference_point) - 5;
                max_reference = min(reference_point) + 40 - 1;
                mask_rfdata = zeros(size(focused_rfdata));
                mask_rfdata(min_reference:max_reference,:) = 1;
                focused_rfdata_masked = mask_rfdata .* focused_rfdata;
                for jj = 1:length(target_element)-1
                    %波面形状評価
                    displacement(target_element(jj),kk,ll,mm)  = - finddelay(focused_rfdata(:,target_element(jj+1)),...
                        focused_rfdata_masked(:,target_element(jj)),max_delay(target_element(jj)));
                end
                displacement_ref(target_element(1:end-1),kk,ll,mm) = diff(reference_point(target_element));
                error_wavefront(kk,ll,mm) = sum(abs(displacement(target_element(1:end-1),kk,ll,mm) - ...
                    displacement_ref(target_element(1:end-1),kk,ll,mm)))/(length(target_element)-1);
                %% 出力データの保存
                dst_path = sprintf('H:/result/2018_12_03_DAS_correct/case26/correct/lateral%0.1fmm/true%d_assumption%d',...
                    focal_point(1,9)*1e3,rate_IMCL(1,ll),rate_IMCL(1,kk));
                if ~exist(dst_path, 'dir')
                    mkdir(dst_path);
                end
                displacement_tmp = displacement(:,kk,ll,mm);
                save([dst_path,'\pre_rfdata.mat'],'focused_rfdata','focused_rfdata_masked','displacement_tmp','reference_point');
                figure;
                imagesc(focused_rfdata_masked);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),'red');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                xlim([min(target_element) max(target_element)])
                ylim([min_reference max_reference])
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_masked.fig'])
                exportfig([dst_path,'\pre_rfdata_masked'],'png',[400,400])
                close gcf
                figure;
                imagesc(focused_rfdata);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),2,'red');
                scatter(min(target_element):max(target_element),point_max_in_mask(min(target_element):max(target_element)),2,'blue','filled');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_whole.fig'])
                exportfig([dst_path,'\pre_rfdata_whole'],'png',[400,400])
                close gcf
                figure;
                plot(focal_depth*1e3,focal_signal_total(:,kk,ll,mm));
                hold on
                scatter(focal_depth(ind_max_signal(kk,ll,mm))*1e3,max_signal(kk,ll,mm),'red','filled')
                hold off
                xlabel('焦点深さ[mm]')
                ylabel('信号強度[au]')
                savefig([dst_path,'\pre_focal_signal.fig'])
                exportfig([dst_path,'\pre_focal_signal'],'png',[400,400])
                close gcf
                fprintf('prepareing...lateral number is %d, medium number is %d, and estimation number is %d\n',mm,ll,kk);
            end
        end
        %% 境界位置を特定する．(境界位置数：num_medium)
        [~,ind_ind_max_signal] = max(max_signal(:,ll,mm));
        detected_boundary(:,ll,mm) =  ind_max_signal(ind_ind_max_signal,ll,mm);
        %% 各媒質(num_medium)に対して同一の焦点位置でのRFデータから波形形状等を評価する．
        for kk = 1:num_rate_IMCL
                %% 反射強度プロファイルを求める．
            v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
            for ii = detected_boundary(1,ll,mm)
                target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                %受信用の参照点算出
                for jj = 1:num_echo_receiver
                    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                    delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                    reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                end
                %送信ビームフォーミング（共通）
                focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                %受信ビームフォーカシング
                focal_signal(kk,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
            end
                %% 参照点と波面との形状相関を求める．（相互相関法）
            %RFデータマスキング（均質性評価のため)
            min_reference = min(reference_point) - 5;
            max_reference = min(reference_point) + 40 - 1;
            mask_rfdata = zeros(size(focused_rfdata));
            mask_rfdata(min_reference:max_reference,:) = 1;
            focused_rfdata_masked = mask_rfdata .* focused_rfdata;
            for jj = 1:length(target_element)-1
                %波面形状評価
                displacement(target_element(jj),kk,ll,mm)  = - finddelay(focused_rfdata(:,target_element(jj+1)),...
                    focused_rfdata_masked(:,target_element(jj)),max_delay(target_element(jj)));
            end
            displacement_ref(target_element(1:end-1),kk,ll,mm) = diff(reference_point(target_element));
            error_wavefront(kk,ll,mm) = sum(abs(displacement(target_element(1:end-1),kk,ll,mm) - ...
                displacement_ref(target_element(1:end-1),kk,ll,mm)))/(length(target_element)-1);
            dst_path = sprintf('H:/result/2018_12_03_DAS_correct/case26/correct/lateral%0.1fmm/true%d_assumption%d',...
                focal_point(1,9)*1e3,rate_IMCL(1,ll),rate_IMCL(1,kk));
            if ~exist(dst_path, 'dir')
                mkdir(dst_path);
            end
            displacement_tmp = displacement(:,kk,ll,mm);
            save([dst_path,'\rfdata.mat'],'focused_rfdata','focused_rfdata_masked','displacement_tmp','reference_point');
            figure;
            imagesc(focused_rfdata_masked);
            hold on
            scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),'red');
            hold off
            xlabel('receiver[ch]');
            ylabel('time[sample]');
            axis square;
            axis tight;
            xlim([min(target_element) max(target_element)])
            ylim([min_reference max_reference])
            colormap(bone);
            colorbar;
            caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
            savefig([dst_path,'\rfdata_masked.fig'])
            exportfig([dst_path,'\rfdata_masked'],'png',[400,400])
            close gcf
            figure;
            imagesc(focused_rfdata);
            hold on
            scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),2,'red');
            scatter(min(target_element):max(target_element),point_max_in_mask(min(target_element):max(target_element)),2,'blue','filled');
            hold off
            xlabel('receiver[ch]');
            ylabel('time[sample]');
            axis square;
            axis tight;
            colormap(bone);
            colorbar;
            caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
            savefig([dst_path,'\rfdata_whole.fig'])
            exportfig([dst_path,'\rfdata_whole'],'png',[400,400])
            close gcf
            fprintf('detected... lateral number is %d, medium number is %d, and estimation number is %d\n',mm,ll,kk);
        end
    end
    %% 出力データの保存
    max_signal_normalized = zeros(num_rate_IMCL,num_medium,num_lateral);
    for ll = 1:num_medium
        tmp2 = focal_signal(:,ll,mm);
        max_signal_normalized(:,ll,mm) = tmp2/max(tmp2);
    end
    
    dst_path = sprintf('H:/result/2018_12_03_DAS_correct/case26/correct/lateral%0.1fmm/',focal_point(1,9)*1e3);
    save([dst_path,'\result.mat'],'ind_max_signal','focal_depth','displacement',...
        'reference_point','focal_signal_total','focal_signal','max_signal','max_signal_normalized');
    
    
    figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,ind_max_signal(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('boundary depth[mm]')
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\pre_boundary_eng'],'png',[300,250]);
    close gcf
    
    figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,ind_max_signal(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('媒質境界深度[mm]')
    xlabel('正解IMCL占有率 [%]');
    ylabel('予測IMCL占有率 [%]');
    exportfig([dst_path,'\pre_boundary'],'png',[300,250]);
    close gcf
    
        figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,detected_boundary(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('boundary depth[mm]')
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\detected_boundary_eng'],'png',[300,250]);
    close gcf
    
    figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,detected_boundary(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('媒質境界深度[mm]')
    xlabel('正解IMCL占有率 [%]');
    ylabel('予測IMCL占有率 [%]');
    exportfig([dst_path,'\detected_boundary'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(error_wavefront(:,:,mm));
    title('wavefront distortion');
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]')
    exportfig([dst_path,'\detected_wavefront_eng'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(error_wavefront(:,:,mm));
    title('波面歪み指標');
    xlabel('正解IMCL占有率 [%]');
    ylabel('予測IMCL占有率 [%]');
    exportfig([dst_path,'\detected_wavefront]]'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(max_signal(:,:,mm));
    title('信号強度[au]');
    xlabel('正解IMCL占有率 [%]');
    ylabel('予測IMCL占有率 [%]');
    exportfig([dst_path,'\pre_intensity'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(max_signal(:,:,mm));
    title('Intensity[au]');
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\pre_intensity_eng'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(focal_signal(:,:,mm));
    title('信号強度[au]');
    xlabel('正解IMCL占有率 [%]');
    ylabel('予測IMCL占有率 [%]');
    exportfig([dst_path,'\detected_intensity'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(focal_signal(:,:,mm));
    title('Intensity[au]');
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\detected_intensity_eng'],'png',[300,250]);
    close gcf    
end
%% すべての変数を保存
savefilename = sprintf('H:/result/2018_12_03_DAS_correct/case26/correct/total_result.mat');
save(savefilename);