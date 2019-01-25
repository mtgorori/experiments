load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/kgrid.mat")
load("H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL1.0_pure/rfdata.mat")
casename = [1,2,17,26,52];
% 音速値
v_fat        = 1450;%[m/s]
v_muscle = 1580;%[m/s]
den_muscle = 1040;%筋肉の密度[kg/m3]
den_fat = 920;%脂肪の密度[kg/m3]
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

load("H:\result\2018_12_19_IMCL_direct_estimation_multi_element\case26\2018_12_25_variousF&lateral\2018_12_28\case26IMCL0%\result.mat")

wire_point = assumed_point(:,ind_estimate_d,ind_estimate_l);
ind_wire_point = zeros(1,2);
[~,ind_wire_point(1)] = min(abs(kgrid.x_vec - wire_point(2)));
[~,ind_wire_point(2)] = min(abs(kgrid.y_vec - wire_point(1)));
load("H:/data/kwave/medium/2018_09_28_realisticScatter_variousIMCL/corrected/case26_IMCL0.0_pure.mat")
size_medium = size(medium.sound_speed);
medium.sound_speed = ones(size_medium)*v_muscle;
medium.sound_speed(ind_wire_point(2),ind_wire_point(1)) = v_fat;
medium.density = ones(size_medium)*den_muscle;
medium.density(ind_wire_point(2),ind_wire_point(1)) = den_fat;

save("H:\data\kwave\medium\2019_01_23_wire_case26_maxcorr\IMCL0%",'medium','kgrid','param')
