%%%%%%%%%%%%%%%%%%%%
% 対象：case26
% 波面遅延プロファイルにより音速推定
% 仮定遅延プロファイルと実測遅延プロファイルの相互相関
% を用いて音速を推定する
% 媒質   IMCL0%
%%%%%%%%%%%%%%%%%%%%

%% 初期設定（共通）
clear
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/medium/2018_09_28_realisticScatter_variousIMCL/corrected/case26_IMCL0.0_pure.mat")
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/rfdata.mat")
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/kgrid.mat")
load("H:/experiments/2018_12_19_IMCL_direct_estimation_multi_element/case26/condition/F&lateralchange/2018_12_28_case26.mat")

% 音速値
v_fat        = 1450;%[m/s]
v_muscle = 1580;%[m/s]

% IMCL割合（正解）
IMCL_rate                  = linspace(0,20,11);%[%]
num_IMCL                  = length(IMCL_rate);
v_muscle_with_IMCL = v_fat * IMCL_rate/100 + v_muscle*(1-IMCL_rate/100);%正解音速[m/s]

% 素子配置
t_facing_distance      = 0.04;%[m]
[~,num_receiver,~]  = size(rfdata);
num_transmitter        = num_receiver/2;
num_receiver             = num_receiver/2;
element_pitch           = abs(t_pos(1,1) - t_pos(1,2));
minimum_elementNum = 20;
lateral_range_max = 4.0*1e-3;
lateral_range_min = -4.0*1e-3;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mmまで
lateral_focus_point = linspace(lateral_range_min,lateral_range_max,num_lateral);

% 探索位置
assumed_depth          = 19e-3:-kgrid.dx:0;
assumed_distance      = 20e-3 - assumed_depth;
num_assumed_depth = length(assumed_depth);
ind_assumed_depth   = zeros(num_assumed_depth,1);% kgrid上では境界位置はどのインデックスで表されるかをもとめる．
assumed_point         = zeros(2,num_assumed_depth,num_lateral);
for i = 1:num_assumed_depth
    ind_assumed_depth(i) = find(single(kgrid.x_vec) == single(assumed_depth(i)));
end
for i = 1:num_lateral
    assumed_point(1,:,i) = lateral_focus_point(i);% 送信フォーカス点
    assumed_point(2,:,i) = kgrid.x_vec(ind_assumed_depth);
end

% 遅延プロファイル
distance_from_assumed_point = zeros(1,num_transmitter);
distance_round_trip                  = zeros(1,num_transmitter);
delay_time_assumed                = zeros(1,num_transmitter);

% 推定値
assumed_IMCL_rate                  = linspace(1,20,20);%[%]
num_assumed_IMCL                  = length(assumed_IMCL_rate);
v_muscle_with_assumed_IMCL = v_fat * assumed_IMCL_rate/100 + v_muscle*(1-assumed_IMCL_rate/100);%正解音速[m/s]
assumed_SOS = linspace(v_muscle,v_muscle_with_assumed_IMCL(end),num_assumed_IMCL);
num_assumed_SOS  = length(assumed_SOS);
correlation                = 0;
estimated_velocity   = zeros(1,num_IMCL);
estimated_IMCL       = zeros(1,num_IMCL);
best_lateral              = zeros(1,num_IMCL);

% 正解値
correct_velocity = zeros(1,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end
% for nn = 1:num_IMCL-1
%     loadfilename = sprintf("H:/result/2018_12_19_IMCL_direct_estimation_multi_element/case26/2018_12_25_variousF&lateral/2018_12_28/case26IMCL%d%%/result.mat",IMCL_rate(nn));
%     load(loadfilename);
%     estimated_velocity(1,nn) = assumed_SOS(1,ind_estimate_v);
%     estimated_IMCL(1,nn) = 100*((v_muscle-estimated_velocity(1,nn))/(v_muscle-v_fat));
%     best_lateral(1,nn) = lateral_focus_point(1,ind_estimate_l);
% end
%% 音速推定処理部
% 仮定遅延プロファイルと実測遅延プロファイルの相互相関を求める
for nn = 1
    load("H:\result\2018_12_19_IMCL_direct_estimation_multi_element\case26\2018_12_25_variousF&lateral\2018_12_28\case26IMCL0%\result.mat")
    correlation = 0;
    loadpath = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',IMCL_rate(nn));
    load([loadpath,'/rfdata.mat'])
    load([loadpath,'/kgrid.mat'])
    load([loadpath,'/sourse_wave.mat'])
    [num_sample,~,~] = size(rfdata);
    
    % rfデータ整形%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % サンプル数を4倍にするためにスプライン補間
    interp_sourcewave                     = interp1(source_wave,linspace(1,num_sample,num_sample*4)','spline');
    [~,offset_interp_sourcewave]   = max(interp_sourcewave);% 送信波形の最大値を取る点．遅延曲線に散布図を重ね合わせることに使う．
    rfdata_echo_only                       = zeros(num_sample,num_receiver,num_transmitter);
    
    for ii = 1:num_transmitter
        for jj = 1:num_receiver
            delay_transmitted_wave           = round(((abs(t_pos(1,jj) - t_pos(1,ii)))/v_muscle)/(kgrid.dt));
            rfdata_echo_only(:,jj,ii)            = rfdata(:,jj,ii);
            rfdata_echo_only(1:delay_transmitted_wave+50,jj,ii) = 0;
        end
    end
    
    % 音速推定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for mm = ind_estimate_l
        for  kk = ind_estimate_d
            % 駆動素子の選択
            target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,kk,mm);
            
            for ll = ind_estimate_v
                
                % 送信フォーカス
                [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,mm,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver);
                load("H:\data\kwave\result\2019_01_23_wire_case26_maxcorr\IMCL0%\rfdata.mat");
                for ii = 1:num_transmitter
                    for jj = 1:num_receiver
                        delay_transmitted_wave           = round(((abs(t_pos(1,jj) - t_pos(1,ii)))/v_muscle)/(kgrid.dt));
                        rfdata_echo_only(:,jj,ii)            = rfdata(:,jj,ii);
                        rfdata_echo_only(1:delay_transmitted_wave+50,jj,ii) = 0;
                    end
                end
                [reference_rfdata,~] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,mm,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver);
                interp_t_array = interp1(kgrid.t_array,linspace(1,num_sample,num_sample*4),'spline');
                distance_round_trip   = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
                delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
                [~,ind_cutoff] = min(abs(5/1e6 - interp_t_array));
                for iii = 1:num_receiver
                    offsets = delay_time_assumed(iii)-min(delay_time_assumed);
                    reference_rfdata(ind_cutoff+offsets:end,iii) = 0;
                end
                % 自己相関係数の最大値を各信号で求める（rfdata, source wave）
                auto_correlation_rfdata            = diag(focused_rfdata.' * focused_rfdata);
                auto_correlation_reference = diag(reference_rfdata.' * reference_rfdata);
                
                %                 % 遅延プロファイルの仮定
                %                 distance_round_trip   = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
                %                 delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
                %                 delay_sourcewave     = zeros(num_sample*4, num_receiver);
                
                for ii = 1:length(target_element)
                    % 相互相関の積算
                    correlation = correlation+ (focused_rfdata(:,target_element(ii)).'*reference_rfdata(:,target_element(ii))...
                        /sqrt(auto_correlation_rfdata(target_element(ii),1)*auto_correlation_reference(target_element(ii),1)));
                end
                
                correlation = correlation/length(target_element);
                
            end
            
            %             dispname = sprintf('estimated depth # is %d, lateral # is %0.1f, IMCL # is %d, target element # is %d',kk,mm,nn,length(target_element));
            %             disp(dispname) %#ok<DSPS>
            
        end
    end
    
    figure;
    plot(interp_t_array*1e6,focused_rfdata(:,47)/max(abs(focused_rfdata(:,47))));
    hold on
    plot(interp_t_array*1e6,reference_rfdata(:,47)/max(abs(reference_rfdata(:,47))));
    hold off
    
    xlabel('時刻[μs]')
    ylabel('音圧[a.u.]')
    legend('実データ','参照データ')
    exportfig("H:\result\2019_01_23_IMCL_referencial_wave\comp_reference_from_wiredata",'png',[300,200])
    xlim([3 5])
    exportfig("H:\result\2019_01_23_IMCL_referencial_wave\comp_reference_from_wiredata_detail",'png',[300,200])
    
    interp_sourcewave                     = interp1(source_wave,linspace(1,num_sample,num_sample*4)','spline');
    [~,offset_interp_sourcewave]   = max(interp_sourcewave);% 送信波形の最大値を取る点．遅延曲線に散布図を重ね合わせることに使う．
    
    
    delay_sourcewave     = zeros(num_sample*4, num_receiver);
    
    for ii = 1:length(target_element)
        source_wave2cat = interp_sourcewave;
        source_wave2cat(num_sample*4-delay_time_assumed(target_element(ii))+1:end,1)=NaN;
        source_wave2cat(isnan(source_wave2cat)) = [];
        delay_sourcewave(:,target_element(ii)) = cat(1,zeros(delay_time_assumed(target_element(ii)),1),source_wave2cat);
    end
    figure;
    plot(interp_t_array*1e6,focused_rfdata(:,47)/max(abs(focused_rfdata(:,47))));
    hold on
    plot(interp_t_array*1e6,delay_sourcewave(:,47)/max(abs(delay_sourcewave(:,47))));
    hold off
    xlabel('時刻[μs]')
    ylabel('音圧[a.u.]')
    legend('実データ','参照データ')
    exportfig("H:\result\2019_01_23_IMCL_referencial_wave\comp_reference_from_sourcewave",'png',[300,200])
    xlim([3 5])
    exportfig("H:\result\2019_01_23_IMCL_referencial_wave\comp_reference_from_sourcewave_detail",'png',[300,200])
end


%% 関数部
function target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,kk,mm)
target_element = find((-assumed_distance(1,kk) + assumed_point(1,1,mm)<=t_pos(1,1:num_receiver))&...
    ((t_pos(1,1:num_receiver)<=assumed_distance(1,kk)+ assumed_point(1,1,mm))));
if length(target_element) < minimum_elementNum
    [~,ind_central_element] = min(abs(t_pos(1,:)-lateral_focus_point(mm)));
    target_element = ceil(ind_central_element-minimum_elementNum/2+1):ceil(ind_central_element+minimum_elementNum/2);
end
end


function [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,mm,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver)
distance_from_assumed_point = zeros(1,num_transmitter);

for ii = 1:num_transmitter
    distance_from_assumed_point(1,ii) = norm(assumed_point(:,kk,mm)-t_pos(:,ii));%[m]
end

delay_length_focusing = distance_from_assumed_point - min(distance_from_assumed_point);
delay_time_focusing    = round(delay_length_focusing/ assumed_SOS(ll) / (kgrid.dt));
focused_rfdata        = transmit_focus(rfdata_echo_only,target_element,...
    num_sample,delay_time_focusing,num_receiver);
interp_rfdata = zeros(num_sample*4,num_receiver);

for jj = 1:num_receiver
    interp_rfdata(:,jj) = interp1(focused_rfdata(:,jj),linspace(1,num_sample,num_sample*4),'spline');
end

focused_rfdata = interp_rfdata;
end