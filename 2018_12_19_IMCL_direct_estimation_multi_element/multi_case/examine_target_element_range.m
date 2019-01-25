%% 初期設定（共通）
clear
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/medium/2018_09_28_realisticScatter_variousIMCL/corrected/case26_IMCL0.0_pure.mat")
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/rfdata.mat")
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/kgrid.mat")
case_name = [1,2,17,26,52];
num_case = length(case_name);
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
correlation                = zeros(num_assumed_SOS,num_assumed_depth,num_lateral);
estimated_velocity   = zeros(1,num_IMCL);
estimated_IMCL       = zeros(1,num_IMCL);
best_lateral              = zeros(1,num_IMCL);

% 正解値
correct_velocity = zeros(1,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end

max_target_element = zeros(num_assumed_depth,num_lateral);
min_target_element = zeros(num_assumed_depth,num_lateral);
max_depth = zeros(num_case,num_IMCL);
%% 音速推定処理部
% 仮定遅延プロファイルと実測遅延プロファイルの相互相関を求め
for pp = 1:num_case
    for nn = 1:num_IMCL
        loadfilename = sprintf('H:/result/2018_12_19_IMCL_direct_estimation_multi_element/multi_case/2018_12_28_variousF&lateral/case%dIMCL%d%%/result.mat',case_name(pp),IMCL_rate(nn));
        load(loadfilename,'ind_estimate_d');
        max_depth(pp,nn) = ind_estimate_d;
    end
end
limit_d = max(max(max_depth));
for kk = 1:limit_d+5
    for mm = 1:num_lateral
        target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,kk,mm);
        max_target_element(kk,mm) = max(target_element);
        min_target_element(kk,mm) = min(target_element);
    end
end
max_element = max(max(max_target_element));
min_element = min(min(min_target_element));

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