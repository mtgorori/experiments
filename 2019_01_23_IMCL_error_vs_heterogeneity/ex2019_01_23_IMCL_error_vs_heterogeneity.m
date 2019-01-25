%%%%%%%%%%%%%%%%%%
%各ケースでの推定誤差平均値と平均不均質度をまとめてプロット．
%相関を求める．
%%%%%%%%%%%%%%%%%%

%% 初期設定
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/kgrid.mat")
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/rfdata.mat")
casename = [1,2,17,26,52];
% 音速値
v_fat        = 1450;%[m/s]
v_muscle = 1580;%[m/s]

% IMCL割合（正解）
IMCL_rate                  = linspace(0,20,11);%[%]
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
num_case = length(casename);
num_IMCL = length(IMCL_rate);

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
correlation                = zeros(num_assumed_SOS,num_assumed_depth,num_lateral);

best_lateral              = zeros(1,num_IMCL);

% 正解値
correct_velocity = zeros(1,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end


fat_in_region = zeros(num_case,num_IMCL);
estimation_error = zeros(num_case,num_IMCL);
best_ind = zeros(2,1) ;
start_ind = zeros(2,1);
end_ind = zeros(2,1);
for ii = 1:num_case
    
    for jj = 1:num_IMCL
        if ii == 4
            loadpath = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',IMCL_rate(jj));
            load([loadpath,'/medium.mat'])
        else
            loadpath = sprintf('H:/data/kwave/result/2018_12_28_other_case_variousIMCL/case%d_IMCL%0.1f',casename(ii),IMCL_rate(jj));
            load([loadpath,'/medium.mat'])
        end
        loadfilename = sprintf("H:/result/2018_12_19_IMCL_direct_estimation_multi_element/multi_case/2018_12_28_variousF&lateral/case%dIMCL%d%%/result.mat",casename(ii),IMCL_rate(jj));
        load(loadfilename);
        target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,ind_estimate_d,ind_estimate_l);
        start_element = target_element(1);
        end_element = target_element(end);
        start_element_pos = t_pos(:,start_element);
        end_element_pos = t_pos(:,end_element);
        best_point = single(assumed_point(:,ind_estimate_d,ind_estimate_l));
        [~,start_ind(2)] = min(abs(kgrid.y_vec  -  start_element_pos(1)));%全部indexに置き換えたい
        [~,end_ind(2)] = min(abs(kgrid.y_vec -  end_element_pos(1)));
        start_ind(1) = find(kgrid.y_vec == t_pos(2,1));
        end_ind(1) = find(kgrid.y_vec == t_pos(2,1));
        [~,best_ind(1)] = min(abs(kgrid.y_vec - best_point(2)));
        [~,best_ind(2)] = min(abs(kgrid.y_vec - best_point(1)));
        lower_line = zeros(2,start_ind(1)-best_ind(1)+1);
        lower_line(1,:) = linspace(best_ind(1),start_ind(1),start_ind(1)-best_ind(1)+1);
        lower_line(2,:) = round(linspace(best_ind(2),start_ind(2),start_ind(1)-best_ind(1)+1));
        upper_line = zeros(2,end_ind(1)-best_ind(1)+1);
        upper_line(1,:) = linspace(best_ind(1),end_ind(1),end_ind(1)-best_ind(1)+1);
        upper_line(2,:) = round(linspace(best_ind(2),end_ind(2),end_ind(1)-best_ind(1)+1));
        total_pixel = 0;
        fat_pixel = 0;
        for kk = 1:start_ind(1)-best_ind(1)+1
            total_pixel = total_pixel + (upper_line(2,kk)-lower_line(2,kk)+1);
            fat_pixel = fat_pixel + sum(medium.sound_speed(upper_line(1,kk),lower_line(2,kk):upper_line(2,kk)) == v_fat);
        end
        fat_in_region(ii,jj) = fat_pixel / total_pixel * 100;
        estimated_velocity = assumed_SOS(1,ind_estimate_v);
        estimated_IMCL = 100*((v_muscle-estimated_velocity)/(v_muscle-v_fat));
        estimation_error(ii,jj) = abs(estimated_IMCL-IMCL_rate(jj));
    end
end

%% 推定誤差
ave_estimation_error = zeros(1,num_case);
ave_fat_in_region = zeros(1,num_case);

for ii = 1:num_case
    ave_estimation_error(1,ii) = mean(estimation_error(ii,:));
    ave_fat_in_region(1,ii) = mean(fat_in_region(ii,:));
end

%% 図示
figure;
scatter(ave_fat_in_region,ave_estimation_error);
xlabel('伝搬領域内でのEMCL含有率の平均[%]');
ylabel('IMCL含有率推定誤差の平均[%]');
mdl = fitlm(ave_fat_in_region,ave_estimation_error);
axis square
exportfig("H:\result\2019_01_23_IMCL_error_vs_heterogeneity\EMCL_vs_error",'png',[200,200]);

%% 関数
function target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,kk,mm)
target_element = find((-assumed_distance(1,kk) + assumed_point(1,1,mm)<=t_pos(1,1:num_receiver))&...
    ((t_pos(1,1:num_receiver)<=assumed_distance(1,kk)+ assumed_point(1,1,mm))));
if length(target_element) < minimum_elementNum
    [~,ind_central_element] = min(abs(t_pos(1,:)-lateral_focus_point(mm)));
    target_element = ceil(ind_central_element-minimum_elementNum/2+1):ceil(ind_central_element+minimum_elementNum/2);
end
end